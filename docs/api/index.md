# popnei

popnei is a library for population genetics calculations.

It reads and writes VCF files and its own variants file (.nei).
With popnei you can:

- filter variants or individuals, for instance by missing data or MAF;
- calculate basic statistics, such as observed or expected heterozygosity,
  the diversity of each population or F_IS;
- calculate distances between individuals and between populations;
- calculate the kinship between individuals;
- do principal component and principal coordinate analyses;
- analyze linkage disequilibrium; and
- run GWAS.

popnei does not load the whole variant file into memory, so it can
process an arbitrary number of variants.

This site documents popnei's Python package.

The [getting started](getting_started.md) tutorial goes through the
installation and a standard population genetics analysis.

```{toctree}
:maxdepth: 1
:hidden:

getting_started
```

```{toctree}
:caption: Reference
:maxdepth: 1

reference/variants
reference/stats
reference/diversity
reference/dists
reference/pca
reference/kinship
reference/ld
reference/gwas
```
