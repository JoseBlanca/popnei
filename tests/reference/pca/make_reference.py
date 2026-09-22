"""The reference data of the PCA, docs/specs/pca.md.

Run from a checkout of pyNei at the commit that pyproject.toml names, with
this directory as the argument:

    PYTHONPATH=. uv run python <popnei>/tests/reference/pca/make_reference.py <popnei>/tests/reference/pca
    Rscript <popnei>/tests/reference/pca/reference.R <popnei>/tests/reference/pca

It writes the panel of pyNei's test/gwas_reference/sim_missing.vars as a
VCF, which popnei reads, the dosage matrices that R reads, the iris table,
and what pyNei gives for each. reference.R then writes what R's prcomp
gives. The sign of every component is set as the spec says: the projection
with the largest absolute value is positive.
"""

import sys
from pathlib import Path

import numpy
import pandas

from pynei import load_vars
from pynei.pca import create_012_gt_matrix, do_pca, do_pca_from_variants
from pynei.variants import Variants
from test.datasets import IRIS

OUT = Path(sys.argv[1])

WORKED = {
    "v0": "0/0 0/1 1/1 0/0 0/1",
    "v1": "1/1 1/1 0/1 ./. 1/1",
    "v2": "0/0 0/0 0/0 0/0 0/0",
    "v3": "0/1 0/1 0/1 0/1 0/1",
    "v4": "0/2 0/0 2/2 0/. 0/0",
    "v5": "0/1 1/2 0/0 0/0 2/2",
}


def fix_signs(projections, princomps=None):
    signs = numpy.sign(
        projections.values[numpy.abs(projections.values).argmax(axis=0), range(projections.shape[1])]
    )
    projections = projections * signs
    if princomps is not None:
        princomps = princomps.mul(signs, axis=0)
    return projections, princomps


def write(name, res, num_comps):
    projections, princomps = fix_signs(res.projections, res.princomps)
    projections.iloc[:, :num_comps].to_csv(OUT / f"{name}.pynei.projections.tsv", sep="\t", float_format="%.12g")
    res.explained_variance_percent.iloc[:num_comps].to_csv(OUT / f"{name}.pynei.percent.tsv", sep="\t", float_format="%.12g", header=False)
    princomps.iloc[:num_comps].to_csv(OUT / f"{name}.pynei.princomps.tsv", sep="\t", float_format="%.12g")


def gts_of(lines):
    rows = [[[-1 if a == "." else int(a) for a in gt.split("/")] for gt in line.split()] for line in lines]
    return numpy.array(rows)


def write_vcf(path, gts, individuals):
    with open(path, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n##contig=<ID=1>\n")
        vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(individuals) + "\n")
        for idx, var in enumerate(gts):
            num_alt = max(int(var.max()), 1)
            alts = ",".join("CGT"[:num_alt])
            calls = "\t".join("/".join("." if a < 0 else str(a) for a in gt) for gt in var)
            vcf.write(f"1\t{idx + 1}\tvar{idx:04d}\tA\t{alts}\t.\t.\t.\tGT\t{calls}\n")


def mat012(variants, **kwargs):
    return create_012_gt_matrix(variants, **kwargs)


# the panel
panel = load_vars("test/gwas_reference/sim_missing.vars")
individuals = list(panel.samples)
gts = numpy.vstack([chunk.gts.gt_values for chunk in load_vars("test/gwas_reference/sim_missing.vars").iter_vars_chunks()])
write_vcf(OUT / "sim_missing.vcf", gts, individuals)
numpy.savetxt(OUT / "sim_missing.mat012.tsv", mat012(load_vars("test/gwas_reference/sim_missing.vars")), fmt="%d", delimiter="\t")
write("sim_missing", do_pca_from_variants(load_vars("test/gwas_reference/sim_missing.vars")), 10)

# the worked example, without its variant of three alleles and with it.
# pyNei counts the alleles of the whole chunk, so it refuses the first one
# too, where no variant has more than two alleles but the chunk holds 0, 1
# and 2. transform_to_biallelic=True does not change the dosage of a variant
# of two alleles, so it is passed in both.
names = ["i0", "i1", "i2", "i3", "i4"]
for name, keys, kwargs in (("worked", ["v0", "v1", "v2", "v3", "v4"], {"transform_to_biallelic": True}), ("worked3", list(WORKED), {"transform_to_biallelic": True})):
    gts = gts_of([WORKED[key] for key in keys])
    write_vcf(OUT / f"{name}.vcf", gts, names)
    numpy.savetxt(OUT / f"{name}.mat012.tsv", mat012(Variants.from_gt_array(gts, samples=names), **kwargs), fmt="%d", delimiter="\t")
    write(name, do_pca_from_variants(Variants.from_gt_array(gts, samples=names), **kwargs), 4)

# iris, for the PCA of a table
iris = IRIS["characterization"]
iris.to_csv(OUT / "iris.tsv", sep="\t")
write("iris", do_pca(iris), 4)
write("iris_not_standardized", do_pca(iris, standarize_data=False), 4)
