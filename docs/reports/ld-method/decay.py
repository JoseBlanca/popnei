"""The pairs of each population grouped by their exact distance, for decay.R."""
import os
import subprocess

import numpy

SP = os.environ.get("LD_WORK", ".")
VCF = SP + "/ld.vcf"
MIN_DIST, MAX_DIST = 1, 250000
raw = []; chroms = []; poss = []
for line in open(VCF):
    if line.startswith("#"):
        if line.startswith("#CHROM"): inds = line.rstrip("\n").split("\t")[9:]
        continue
    f = line.rstrip("\n").split("\t")
    chroms.append(f[0]); poss.append(int(f[1]))
    raw.append([[None if a == "." else int(a) for a in c.replace("|", "/").split("/")]
                for c in f[9:]])
poss = numpy.array(poss); NV = len(raw)

def maf_of(v, idxs):
    c = {}
    for i in idxs:
        for a in raw[v][i]:
            if a is not None: c[a] = c.get(a, 0) + 1
    t = sum(c.values()); return max(c.values()) / t if t else None

def run(tag, idxs, max_maf):
    keep = [v for v in range(NV) if (m := maf_of(v, idxs)) is not None and m <= max_maf]
    with open(f"{SP}/{tag}.keep", "w") as fh:
        for i in idxs: fh.write(f"{inds[i]}\t{inds[i]}\n")
    with open(f"{SP}/{tag}.extract", "w") as fh:
        for v in keep: fh.write(f"v{v:04d}\n")
    subprocess.run(["plink2", "--vcf", VCF, "--double-id", "--allow-extra-chr",
                    "--keep", f"{SP}/{tag}.keep", "--extract", f"{SP}/{tag}.extract",
                    "--r2-unphased", "square", "bin", "--out", f"{SP}/{tag}"],
                   capture_output=True, check=True)
    M = numpy.fromfile(f"{SP}/{tag}.unphased.vcor2.bin",
                       dtype=numpy.float64).reshape(len(keep), len(keep))
    acc = {}
    for a in range(len(keep)):
        for b in range(a + 1, len(keep)):
            i, j = keep[a], keep[b]
            if chroms[i] != chroms[j]: continue
            d = int(poss[j] - poss[i])
            if d < MIN_DIST or d > MAX_DIST: continue
            v = float(M[a, b])
            if v != v: continue
            n, s = acc.get(d, (0, 0.0)); acc[d] = (n + 1, s + v)
    with open(f"{SP}/{tag}.decay.tsv", "w") as fh:
        fh.write("dist\tnum_pairs\tsum_r2\n")
        for d in sorted(acc):
            n, s = acc[d]; fh.write(f"{d}\t{n}\t{s!r}\n")
    print(f"{tag}: {len(idxs)} individuals, {len(keep)} of {NV} variants with "
          f"maf<={max_maf}, {len(acc)} distances holding "
          f"{sum(n for n, _ in acc.values())} pairs")

run("all", list(range(100)), 0.95)
run("pop_a", list(range(0, 50)), 0.8)
run("pop_b", list(range(50, 100)), 0.8)
