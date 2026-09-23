"""The score test p-values of every variant, from each candidate null fit,
against GMMAT's glmm.score on pyNei's reference panel."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy, pandas
from fits import fit_a, fit_b, fit_c
import pynei.gwas as g
from pynei import load_vars, calc_kinship

REF = pathlib.Path("/Users/jose/devel/pynei/test/gwas_reference")
variants = load_vars(REF / "sim.vars")
pheno = pandas.read_csv(REF / "phenotypes.csv", index_col="IID")
kin = calc_kinship(variants).matrix.to_numpy()
y = pheno["binom"].to_numpy(float)
design = numpy.column_stack([numpy.ones(len(y)), pheno[["cov1", "cov2"]].to_numpy(float)])
ref = pandas.read_csv(REF / "gmmat_glmm_score.tsv", sep="\t")

# the dosages of every variant, in the order the file has them
chunks = [g._calc_dosages(chunk, None) for chunk in variants.iter_vars_chunks()]
dosages = numpy.concatenate([d for d, _, _ in chunks])
is_poly = numpy.concatenate([p for _, _, p in chunks])
print("variants", dosages.shape[0], " varying", int(is_poly.sum()))

for name, fit in (("A pyNei", fit_a), ("B cholesky", fit_b), ("C identity", fit_c)):
    res = fit(y, design, kin.copy())
    x = dosages[is_poly]
    num = x @ res["py"]
    den = numpy.einsum("ij,ij->i", x @ res["projection"], x)
    p_value = g._chi2_sf_1df(num * num / den)
    variance = 1 / (1 / den)          # the variance of the score, GMMAT's VAR
    p_all = numpy.full(dosages.shape[0], numpy.nan)
    p_all[is_poly] = p_value
    var_all = numpy.full(dosages.shape[0], numpy.nan)
    var_all[is_poly] = den
    log10 = numpy.nanmax(numpy.abs(numpy.log10(p_all / ref["PVAL"].to_numpy())))
    var_rel = numpy.nanmax(numpy.abs(var_all / ref["VAR"].to_numpy() - 1))
    print(f"{name:>11}: tau={res['tau']:.9f}  max |log10(p/p_gmmat)| {log10:.3e}  "
          f"max rel diff of the score variance {var_rel:.3e}")
