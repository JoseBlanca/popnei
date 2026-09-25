import os
import numpy, subprocess
SP = os.environ.get("LD_WORK", ".")
def load(vcf, binfile, n):
    chroms=[];poss=[];raw=[]
    for line in open(vcf):
        if line.startswith("#"): continue
        f=line.rstrip("\n").split("\t")
        chroms.append(f[0]);poss.append(int(f[1]))
        raw.append([[None if a=="." else int(a) for a in c.replace("|","/").split("/")] for c in f[9:]])
    M=numpy.fromfile(binfile,dtype=numpy.float64).reshape(n,n)
    return chroms,numpy.array(poss),raw,M
def has_variance(g):
    ds=[]
    c={}
    for gt in g:
        for a in gt:
            if a is not None: c[a]=c.get(a,0)+1
    if not c: return False
    top=max(c.values()); major=min(a for a,k in c.items() if k==top)
    for gt in g:
        if all(a is not None for a in gt): ds.append(sum(a!=major for a in gt))
    return len(ds)>0 and any(d!=ds[0] for d in ds)
def prune(chroms,poss,raw,M,max_dist,thr):
    kept=[]
    for i in range(len(raw)):
        if not has_variance(raw[i]): continue
        win=[j for j in kept if chroms[j]==chroms[i] and poss[i]-poss[j]<=max_dist]
        if not any(M[j,i]>thr for j in win): kept.append(i)
    return kept

print("== the worked example, 5 variants ==")
c,p,r,M=load(SP+"/example.vcf",SP+"/ex.unphased.vcor2.bin",5)
for md,thr in [(5000,0.5),(5000,0.7),(1000,0.5)]:
    k=prune(c,p,r,M,md,thr)
    print(f"  max_dist {md}, max_allowed_r2 {thr}: kept {[f'v{i+1}' for i in k]}")

print("\n== ld.vcf, 500 variants of 100 individuals ==")
c,p,r,M=load(SP+"/ld.vcf",SP+"/ld_r2.unphased.vcor2.bin",500)
for md,thr in [(10000,0.1),(10000,0.3),(50000,0.3),(250000,0.3)]:
    k=prune(c,p,r,M,md,thr)
    ks=set(k)
    bad=[(a,b) for ia,a in enumerate(k) for b in k[ia+1:]
         if c[a]==c[b] and p[b]-p[a]<=md and M[a,b]>thr]
    unex=[i for i in range(500) if i not in ks and has_variance(r[i])
          and not any(M[j,i]>thr for j in k if j<i and c[j]==c[i] and p[i]-p[j]<=md)]
    subprocess.run(["plink2","--vcf",SP+"/ld.vcf","--double-id","--allow-extra-chr",
                    "--indep-pairwise",f"{md//1000}kb",str(thr),"--out",SP+"/f"],capture_output=True)
    pl=sum(1 for _ in open(SP+"/f.prune.in"))
    first_chr2=[i for i in k if c[i]=="chr2"][:1]
    print(f"  max_dist {md:>6}, max_allowed_r2 {thr}: kept {len(k):>3} of 500 | "
          f"kept pairs above the threshold {len(bad)} | dropped with variance and no kept partner {len(unex)} | "
          f"plink2 --indep-pairwise keeps {pl}")
    print(f"      first five kept {[f'{c[i]}:{p[i]}' for i in k[:5]]}, first kept of chr2 {[f'{c[i]}:{p[i]}' for i in first_chr2]}")
