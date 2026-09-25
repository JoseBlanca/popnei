import os
import numpy, subprocess
SP = os.environ.get("LD_WORK", ".")
VCF=SP+"/ld.vcf"
raw=[];chroms=[];poss=[]
for line in open(VCF):
    if line.startswith("#"):
        if line.startswith("#CHROM"): inds=line.rstrip("\n").split("\t")[9:]
        continue
    f=line.rstrip("\n").split("\t")
    chroms.append(f[0]); poss.append(int(f[1]))
    raw.append([[None if a=="." else int(a) for a in c.replace("|","/").split("/")] for c in f[9:]])
poss=numpy.array(poss); NV=len(raw)
def maf_of(v,idxs):
    c={}
    for i in idxs:
        for a in raw[v][i]:
            if a is not None: c[a]=c.get(a,0)+1
    t=sum(c.values()); return max(c.values())/t if t else None
def bins(M, keep, idxs_pop, min_dist, max_dist, num_bins):
    w=(max_dist-min_dist+1)/num_bins
    n=numpy.zeros(num_bins,dtype=numpy.int64); s=numpy.zeros(num_bins); s2=numpy.zeros(num_bins)
    for a in range(len(keep)):
        for b in range(a+1,len(keep)):
            i,j=keep[a],keep[b]
            if chroms[i]!=chroms[j]: continue
            d=int(poss[j]-poss[i])
            if d<min_dist or d>max_dist: continue
            v=M[a,b]
            if v!=v: continue
            k=min(int((d-min_dist)//w),num_bins-1)
            n[k]+=1; s[k]+=v; s2[k]+=v*v
    return n,s,s2,w
def run(tag, idxs, max_maf, min_dist, max_dist, num_bins):
    keep=[v for v in range(NV) if (m:=maf_of(v,idxs)) is not None and m<=max_maf]
    with open(f"{SP}/x.keep","w") as fh:
        for i in idxs: fh.write(f"{inds[i]}\t{inds[i]}\n")
    with open(f"{SP}/x.extract","w") as fh:
        for v in keep: fh.write(f"v{v:04d}\n")
    subprocess.run(["plink2","--vcf",VCF,"--double-id","--allow-extra-chr","--keep",f"{SP}/x.keep",
                    "--extract",f"{SP}/x.extract","--r2-unphased","square","bin","--out",f"{SP}/x"],
                   capture_output=True)
    M=numpy.fromfile(f"{SP}/x.unphased.vcor2.bin",dtype=numpy.float64).reshape(len(keep),len(keep))
    n,s,s2,w=bins(M,keep,idxs,min_dist,max_dist,num_bins)
    print(f"\n{tag}: {len(idxs)} individuals, {len(keep)} of {NV} variants with maf<={max_maf}, "
          f"bin width {w:g} bp")
    for k in range(num_bins):
        lo=int(min_dist+k*w); hi=int(min_dist+(k+1)*w)-1
        mean=s[k]/n[k]; sd=numpy.sqrt(max(s2[k]/n[k]-mean*mean,0.0))
        print(f"  [{lo:>6}, {hi:>6}]  pairs {n[k]:>5}  mean r2 {mean!r}  sd {sd!r}")
run("one population, every individual", list(range(100)), 0.95, 1, 250000, 10)
run("pop_a", list(range(0,50)), 0.8, 1, 250000, 10)
run("pop_b", list(range(50,100)), 0.8, 1, 250000, 10)
