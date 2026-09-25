import os
import numpy
SP = os.environ.get("LD_WORK", ".")
VCF="/Users/jose/devel/popnei/tests/reference/vcf/many.vcf"
raw=[]
for line in open(VCF):
    if line.startswith("#"): continue
    f=line.rstrip("\n").split("\t")[9:]
    raw.append([[None if a=="." else int(a) for a in c.replace("|","/").split("/")] for c in f])
NV,NI=len(raw),len(raw[0])
def dosages(v, count_half_called_alleles=True):
    c={}
    for g in raw[v]:
        if not count_half_called_alleles and any(a is None for a in g): continue
        for a in g:
            if a is not None: c[a]=c.get(a,0)+1
    if not c: return numpy.full(NI,-1,dtype=numpy.int16)
    top=max(c.values()); major=min(a for a,n in c.items() if n==top)
    out=numpy.empty(NI,dtype=numpy.int16)
    for i,g in enumerate(raw[v]):
        out[i]=-1 if any(a is None for a in g) else sum(a!=major for a in g)
    return out
def r2(x,y):
    m=(x>=0)&(y>=0)
    if m.sum()==0: return numpy.nan
    a=x[m].astype(float); b=y[m].astype(float); a-=a.mean(); b-=b.mean()
    sxx=(a*a).sum(); syy=(b*b).sum()
    return numpy.nan if sxx<=0 or syy<=0 else ((a*b).sum()**2)/(sxx*syy)
M=numpy.fromfile(SP+"/many.unphased.vcor2.bin",dtype=numpy.float64).reshape(NV,NV)
for flag in (True,False):
    D=numpy.array([dosages(v,flag) for v in range(NV)])
    diffs=[]; both_nan=0; only_one=0
    for i in range(0,NV,3):
        for j in range(i+1,NV,7):
            mine=r2(D[i],D[j]); theirs=M[i,j]
            if mine!=mine and theirs!=theirs: both_nan+=1; continue
            if (mine!=mine)!=(theirs!=theirs): only_one+=1; continue
            diffs.append(abs(mine-theirs))
    lbl="half called allele counted" if flag else "half called genotype dropped whole"
    print(f"{lbl}: {len(diffs)} pairs, max |diff| {max(diffs):.3g}, "
          f"both nan {both_nan}, nan on one side {only_one}")
multi=[v for v in range(NV) if len({a for g in raw[v] for a in g if a is not None})>2]
print("multiallelic variants in many.vcf:",len(multi))
