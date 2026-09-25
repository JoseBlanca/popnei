"""The table of "Missing genotypes" of docs/specs/ld.md.

Over the 1.4 million pairs off the diagonal of the panel it prints how far
each of the three rules for a missing genotype lands from plink2's r2: the
individual left out of the pair, the missing genotype given the mean dosage
of its variant, and pyNei's, which leaves it in with a dosage of -1.

Reads $LD_WORK/panel_r2.unphased.vcor2.bin, which plink2 --r2-unphased
writes, and runs from the root of the repository.
"""

import os

import numpy
from pynei.io_vcf import vars_from_vcf
from pynei.ld_calc import _calc_rogers_huff_r2

SP = os.environ.get("LD_WORK", ".")
VCF = "tests/reference/dists/panel.vcf.gz"
NUM_VARS, NUM_INDIS = 1200, 200

plink = numpy.fromfile(
    SP + "/panel_r2.unphased.vcor2.bin", dtype=numpy.float64
).reshape(NUM_VARS, NUM_VARS)
chunks = vars_from_vcf(VCF).iter_vars_chunks()
dosages = numpy.concatenate([c.gts.to_012() for c in chunks], axis=0)
called = dosages >= 0
off_diagonal = ~numpy.eye(NUM_VARS, dtype=bool)


def report(label, r2):
    d = numpy.abs(r2 - plink)[off_diagonal]
    print(
        f"{label}: median {numpy.median(d):.3g}, "
        f"99th pct {numpy.percentile(d, 99):.3g}, largest {d.max():.3g}"
    )


# The individual left out of the pair. Each of the six sums runs over the
# individuals called at both variants, and every one of them is a whole
# number, so the products below are exact in float64.
x = numpy.where(called, dosages, 0).astype(numpy.float64)
both = called.astype(numpy.float64)
n = both @ both.T
sxy = x @ x.T
sx = x @ both.T
sy = sx.T
sxx = (x * x) @ both.T
syy = sxx.T
spread_x = n * sxx - sx * sx
spread_y = n * syy - sy * sy
with numpy.errstate(invalid="ignore", divide="ignore"):
    pairwise = numpy.where(
        (spread_x > 0) & (spread_y > 0),
        (n * sxy - sx * sy) ** 2 / (spread_x * spread_y),
        numpy.nan,
    )
report("the individual is left out of the pair", pairwise)

# The missing genotype takes the mean dosage of its variant.
imputed = dosages.astype(numpy.float64)
for i in range(NUM_VARS):
    missing = ~called[i]
    if missing.any():
        imputed[i][missing] = imputed[i][called[i]].mean()
imputed -= imputed.mean(axis=1, keepdims=True)
ss = numpy.einsum("ij,ij->i", imputed, imputed)
with numpy.errstate(invalid="ignore", divide="ignore"):
    report(
        "the genotype takes the mean dosage of its variant",
        (imputed @ imputed.T) ** 2 / numpy.outer(ss, ss),
    )

# pyNei: the genotype is a dosage of -1. _calc_rogers_huff_r2 returns r and
# not r2, whatever its name says, so it is squared here.
pynei = _calc_rogers_huff_r2(dosages, dosages, check_no_mafs_above=None) ** 2
report("pyNei: the genotype is a dosage of -1", pynei)

theirs = plink[off_diagonal]
print(
    f"plink2 r2 of the panel: median {numpy.median(theirs):.4f}, "
    f"NaN in {numpy.isnan(theirs).sum()} of the {theirs.size} pairs"
)
