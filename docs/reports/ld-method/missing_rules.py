import os
import numpy
from pynei.io_vcf import vars_from_vcf
from pynei.ld_calc import _calc_rogers_huff_r2
SP = os.environ.get("LD_WORK", ".")
p = numpy.fromfile(SP+"/panel_rb.unphased.vcor1.bin", dtype=numpy.float64).reshape(1200,1200)
v = vars_from_vcf("/Users/jose/devel/popnei/tests/reference/dists/panel.vcf.gz")
gts = numpy.concatenate([c.gts.to_012() for c in v.iter_vars_chunks()], axis=0)
py = _calc_rogers_huff_r2(gts, gts, check_no_mafs_above=None)
off = ~numpy.eye(1200, dtype=bool)
d = numpy.abs(py - p)[off]
print("pyNei(-1) vs plink2, |r| difference: median %.4f, 99th pct %.4f, max %.4f" %
      (numpy.median(d), numpy.percentile(d, 99), d.max()))
print("plink |r| median %.4f" % numpy.median(numpy.abs(p[off])))
# mean imputation over the whole matrix
x = gts.astype(float)
for i in range(x.shape[0]):
    m = gts[i] < 0
    if m.any(): x[i][m] = x[i][~m].mean()
x -= x.mean(axis=1, keepdims=True)
ss = numpy.einsum("ij,ij->i", x, x)
mi = (x @ x.T)/numpy.sqrt(numpy.outer(ss, ss))
dm = numpy.abs(mi - p)[off]
print("mean-imputed vs plink2: median %.5f, 99th pct %.5f, max %.5f" %
      (numpy.median(dm), numpy.percentile(dm, 99), dm.max()))
