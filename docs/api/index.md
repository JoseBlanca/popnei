# popnei

popnei is a library of population genetics over the variants of a dataset.
It reads a VCF, or its own vars file, filters the variants and the
individuals, and calculates statistics per variant and per population,
distances between individuals and between populations, principal components
and principal coordinates, the kinship, the linkage disequilibrium and
association studies. Every calculation is written in Rust and runs over the
file one block of variants at a time, so a dataset is never in memory as a
whole. This site documents the Python package, whose names and signatures
follow those of pyNei, the library popnei succeeds.

[Getting started](getting_started.md) installs the package and goes from a
VCF to the main analyses. The reference has one page per subject, generated
from the docstrings of the package. Every name on it is imported from
`popnei` itself: `popnei.do_pca`, not `popnei.pca.do_pca`.

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
