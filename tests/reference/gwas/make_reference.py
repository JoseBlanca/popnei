r"""It writes what plink2, GMMAT, rrBLUP and R's glm give for the association study.

Run from the root of the repository:

    uv run python tests/reference/gwas/make_reference.py

The four programs are run on popnei's own VCFs, the two panels of
`docs/specs/kinship.md`, and not on pyNei's vars files, which popnei cannot
read. pyNei ran the same four on the same genotypes when its own reference was
made, `test/gwas_reference/make_reference.py` in its repository, so at the end
this script compares every number it has produced with what pyNei stored, and
says how far apart they are. They should agree to the digit the programs print:
the genotypes are the same and only the path they arrived by differs.

It needs, and says which one is missing when it cannot find it:

- plink2 v2.0.0-a.7.7 and Rscript with R 4.6.1, GMMAT 1.5.0 and rrBLUP 4.6.3,
  the versions every literal of `docs/specs/gwas.md` was taken with.
- /Users/jose/devel/pynei/test/gwas_reference/, for `phenotypes.csv`,
  `causal_vars.csv` and the stored outputs to compare against. The phenotypes
  are simulated together with the genotypes, so they are copied and not made
  again: making them again would need pyNei's simulation and would give the
  same file.

The two panels, both 200 individuals `s000` to `s199` and 1200 biallelic
diploid variants `var0000` to `var1199` on two chromosomes:

- `panel_called`, every genotype called, `tests/reference/kinship/panel_called.vcf.gz`.
- `panel`, 3 in 100 genotypes missing whole, `tests/reference/dists/panel.vcf.gz`.

The trait is simulated from the genotypes with a heritability of 0.5, five
causal variants of effect 0.6, and two covariates, one continuous and one
binary, with the three subpopulations differing in their mean so that the
structure confounds the trait. `cont` is the continuous trait and `binom` the
0/1 one, the top 40 per cent of a liability.

It writes, beside itself:

- `phenotypes.csv` and `causal_vars.csv`, copied from pyNei.
- `plink2.<panel>.glm.linear.tsv` and `.glm.logistic.hybrid.tsv`: the Wald
  tests of the linear and the logistic model with the two covariates and no
  kinship.
- `gmmat.<panel>.<model>.score.tsv`: GMMAT's score test of the linear and the
  logistic mixed model, with the kinship plink2 wrote for `panel_called`, so
  that the mixed models are tested against a kinship that came from neither
  library.
- `gmmat.null_models.tsv`: the variance components and the covariate effects
  of the two null models GMMAT fitted.
- `rrblup.panel_called.lmm.tsv`: rrBLUP's `GWAS` with the variance components
  held at the null, which is what it calls P3D, and the binary covariate as
  its one fixed effect.
- `r.panel_called.glm.score.tsv`: R's `anova(glm, test = "Rao")`, one logistic
  regression per variant, which is the score test with no kinship.
"""

from pathlib import Path
import gzip
import shutil
import subprocess
import sys
import tempfile

import numpy
import pandas

REF_DIR = Path(__file__).parent
KINSHIP_DIR = REF_DIR.parent / "kinship"
DISTS_DIR = REF_DIR.parent / "dists"
PYNEI_REF = Path("/Users/jose/devel/pynei/test/gwas_reference")

PLINK2_VERSION = "v2.0.0-a.7.7"
R_VERSION = "4.6.1"
GMMAT_VERSION = "1.5.0"
RRBLUP_VERSION = "4.6.3"

PANELS = {
    "panel_called": KINSHIP_DIR / "panel_called.vcf.gz",
    "panel": DISTS_DIR / "panel.vcf.gz",
}

# What pyNei stored for the same genotypes, and what this script writes for it.
AGAINST_PYNEI = {
    "plink2.panel_called.glm.linear.tsv": "plink2.cont.glm.linear.tsv",
    "plink2.panel_called.glm.logistic.hybrid.tsv": "plink2.binom.glm.logistic.hybrid.tsv",
    "gmmat.panel_called.lmm.score.tsv": "gmmat_lmm_score.tsv",
    "gmmat.panel_called.glmm.score.tsv": "gmmat_glmm_score.tsv",
    "gmmat.panel.lmm.score.tsv": "gmmat_lmm_score_missing.tsv",
    "gmmat.panel.glmm.score.tsv": "gmmat_glmm_score_missing.tsv",
    "rrblup.panel_called.lmm.tsv": "rrblup_lmm.tsv",
    "r.panel_called.glm.score.tsv": "r_glm_score.tsv",
    "gmmat.null_models.tsv": "r_null_models.tsv",
}

R_SCRIPT = r"""
suppressMessages({library(GMMAT); library(rrBLUP)})
stopifnot(as.character(packageVersion("GMMAT")) == "%(gmmat)s")
stopifnot(as.character(packageVersion("rrBLUP")) == "%(rrblup)s")

pheno <- read.csv("phenotypes.csv")
rownames(pheno) <- pheno$IID

# the kinship is the one plink2 wrote for the panel with every genotype
# called, so that the mixed models are tested with a kinship that comes from
# neither popnei nor pyNei
K <- as.matrix(read.table("panel_called.rel"))
ids <- read.table("panel_called.rel.id", header = TRUE, comment.char = "")
rownames(K) <- ids[[ncol(ids)]]
colnames(K) <- ids[[ncol(ids)]]

lmm <- glmmkin(cont ~ cov1 + cov2, data = pheno, kins = K, id = "IID",
               family = gaussian(link = "identity"))
glmm <- glmmkin(binom ~ cov1 + cov2, data = pheno, kins = K, id = "IID",
                family = binomial(link = "logit"))
for (panel in c("panel_called", "panel")) {
  glmm.score(lmm, infile = panel,
             outfile = paste0("gmmat.", panel, ".lmm.score.tsv"), MAF.range = c(0, 1))
  glmm.score(glmm, infile = panel,
             outfile = paste0("gmmat.", panel, ".glmm.score.tsv"), MAF.range = c(0, 1))
}
null <- data.frame(
  model = c("lmm", "glmm"),
  tau = c(lmm$theta[2], glmm$theta[2]),
  sigma2 = c(lmm$theta[1], glmm$theta[1]),
  intercept = c(lmm$coefficients[1], glmm$coefficients[1]),
  cov1 = c(lmm$coefficients[2], glmm$coefficients[2]),
  cov2 = c(lmm$coefficients[3], glmm$coefficients[3])
)
write.table(null, "gmmat.null_models.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

dosages <- read.csv("dosages.csv", row.names = 1, check.names = FALSE)
dosages <- as.matrix(dosages)[, pheno$IID]

# rrBLUP takes every fixed effect as a factor, so only the binary covariate
# goes in, and its dosages are counted from -1 to 1
geno <- data.frame(marker = rownames(dosages), chrom = 1,
                   pos = seq_len(nrow(dosages)), dosages - 1, check.names = FALSE)
ph <- data.frame(line = pheno$IID, cont = pheno$cont, cov2 = pheno$cov2)
rr <- GWAS(ph, geno, fixed = "cov2", K = K, n.PC = 0, min.MAF = 0, P3D = TRUE,
           plot = FALSE)
write.table(rr, "rrblup.panel_called.lmm.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)

score <- t(sapply(seq_len(nrow(dosages)), function(idx) {
  g <- dosages[idx, ]
  if (var(g) == 0) return(c(NA, NA))
  fit <- glm(binom ~ cov1 + cov2 + g, data = cbind(pheno, g = g), family = binomial)
  an <- anova(fit, test = "Rao")
  c(an["g", "Rao"], an["g", "Pr(>Chi)"])
}))
write.table(data.frame(id = rownames(dosages), score = score[, 1], p = score[, 2]),
            "r.panel_called.glm.score.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
""" % {"gmmat": GMMAT_VERSION, "rrblup": RRBLUP_VERSION}


def refuse_what_is_missing():
    for tool in ("plink2", "Rscript"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} is not on the PATH")
    printed = subprocess.run(["plink2", "--version"], capture_output=True, text=True,
                             check=True).stdout
    if PLINK2_VERSION not in printed:
        sys.exit(f"this script needs plink2 {PLINK2_VERSION}, and found {printed.strip()}")
    printed = subprocess.run(["Rscript", "-e", "cat(R.version.string)"],
                             capture_output=True, text=True, check=True).stdout
    if R_VERSION not in printed:
        sys.exit(f"this script needs R {R_VERSION}, and found {printed.strip()}")
    for path in list(PANELS.values()) + [PYNEI_REF]:
        if not path.exists():
            sys.exit(f"{path} is not there")


def dosages_of(vcf_gz):
    """The alternative allele count of every genotype, variants x individuals.

    A missing genotype is nan. It is what rrBLUP and R's glm are given, and
    only the panel with every genotype called reaches them.
    """
    ids, rows, individuals = [], [], None
    with gzip.open(vcf_gz, "rt") as fhand:
        for line in fhand:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                individuals = fields[9:]
                continue
            ids.append(fields[2])
            rows.append([
                numpy.nan if gt.startswith(".") else float(gt[0]) + float(gt[2])
                for gt in fields[9:]
            ])
    return pandas.DataFrame(rows, index=ids, columns=individuals)


def main():
    refuse_what_is_missing()
    for name in ("phenotypes.csv", "causal_vars.csv"):
        shutil.copy(PYNEI_REF / name, REF_DIR / name)
    pheno = pandas.read_csv(REF_DIR / "phenotypes.csv")

    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        pheno[["IID", "cont", "binom"]].to_csv(work / "pheno.txt", sep="\t", index=False)
        pheno[["IID", "cov1", "cov2"]].to_csv(work / "covar.txt", sep="\t", index=False)
        shutil.copy(REF_DIR / "phenotypes.csv", work / "phenotypes.csv")
        dosages_of(PANELS["panel_called"]).to_csv(work / "dosages.csv")

        for name, vcf in PANELS.items():
            # --make-bed is what GMMAT reads, --make-rel the kinship R uses
            subprocess.run(
                ["plink2", "--vcf", str(vcf.resolve()), "--pheno", "pheno.txt",
                 "--covar", "covar.txt", "--1", "--make-bed", "--make-rel", "square",
                 "--out", name],
                cwd=work, check=True, capture_output=True)
            subprocess.run(
                ["plink2", "--vcf", str(vcf.resolve()), "--pheno", "pheno.txt",
                 "--covar", "covar.txt", "--1", "--glm", "hide-covar",
                 "--out", f"plink2.{name}"],
                cwd=work, check=True, capture_output=True)
            for kind in ("cont.glm.linear", "binom.glm.logistic.hybrid"):
                written = work / f"plink2.{name}.{kind}"
                if written.exists():
                    trait, _, rest = kind.partition(".")
                    shutil.copy(written, REF_DIR / f"plink2.{name}.{rest}.tsv")
            print(f"plink2 ran on {name}")

        (work / "reference.R").write_text(R_SCRIPT)
        subprocess.run(["Rscript", "reference.R"], cwd=work, check=True)
        print("R, GMMAT and rrBLUP ran")
        for name in AGAINST_PYNEI:
            if name.startswith(("gmmat", "rrblup", "r.")):
                shutil.copy(work / name, REF_DIR / name)

    print("\nagainst what pyNei stored for the same genotypes:")
    worst = 0.0
    for mine, theirs in AGAINST_PYNEI.items():
        ours = pandas.read_csv(REF_DIR / mine, sep="\t")
        pynei = pandas.read_csv(PYNEI_REF / theirs, sep="\t")
        numeric = ours.select_dtypes("number")
        shared = [col for col in numeric.columns if col in pynei.columns]
        diff = max(
            (numpy.nanmax(numpy.abs(numeric[col].to_numpy()
                                    - pynei[col].to_numpy(dtype=float))) or 0.0)
            for col in shared
        ) if shared else float("nan")
        worst = max(worst, diff if diff == diff else 0.0)
        print(f"  {mine:<45} {len(shared):>2} columns, largest difference {diff:.3e}")
    print(f"the largest difference anywhere is {worst:.3e}")


if __name__ == "__main__":
    sys.exit(main())
