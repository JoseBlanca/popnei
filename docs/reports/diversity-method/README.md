# Where the numbers of the diversity spec come from

24 September 2026. The programs behind the tables of
`docs/specs/diversity.md`: the alleles each population of a dataset
called, the private ones among them, the variants that vary in it, the
same three taken down to a common number of called alleles, the folded
site frequency spectrum projected to that number, and F_IS. The spec says
what each number means and what was decided from it; this page says what
each file is and how to run it. Nothing here is part of popnei's build or
of its tests, and the plan that builds the module will take what it needs
into `tests/reference/diversity/make_reference.py`.

pyNei has none of these five calculations, so there is no comparison with
it anywhere here, and every check is against a program outside the
project.

The machine was the owner's Apple M5 Pro, macOS 27.0, with R 4.6.1
holding `vegan` 2.7.6, `adegenet` 2.1.11 and `poppr` 2.9.8, and the Python
environment of popnei, 3.14 with the free threading build, numpy 2.5.3.
`scikit-allel` 1.3.13 and `dadi` 2.4.4 are not in that environment:
`scikit-allel` installs into it and warns that it re-enables the global
interpreter lock, and `dadi` does not build there at all, its `nlopt`
dependency failing, so it was run in a Python 3.12 environment made with
`uv venv --python 3.12`. `hierfstat`, which would have checked the
rarefied allele counts and F_IS in one package, does not install:
`RcppParallel` and then its dependency `gaston` fail to build, as
`docs/specs/dists.md` recorded on 23 September 2026.

Run the Python scripts from this directory. They read the panel of
`docs/specs/stats.md`, `tests/reference/stats/panel.vcf.gz` with
`tests/reference/stats/panel_pops_bcftools.txt`, whose paths are written
into `panel.py` and `check_allel.py`.

- `diversity.py` computes the five quantities as the spec defines them,
  from a VCF and a set of populations. It is what the tables of the spec
  were produced with, and it is the thing the outside programs were
  compared against. It is not popnei and it is not fast: it is the
  arithmetic written out so that it can be read beside the spec.
- `worked.py` runs it on the six variants of five individuals of the
  worked example of `docs/specs/filters.md`, with the two populations of
  the worked example of `docs/specs/stats.md`, at no draw and at draws of
  2 and 4 called alleles. Its output is the worked example of every item
  of the spec, and every one of its numbers was also worked out by hand in
  the spec.
- `panel.py <draw>` runs it on the panel, 1200 biallelic diploid variants
  of 200 individuals in three populations, and writes `panel_counts.tsv`,
  the allele counts of each population at each variant, which the R checks
  and the `dadi` check read. The spec's tables are its output at a draw of
  20.
- `check_vegan.R` gives the rarefied number of alleles of each population
  with `vegan::rarefy`, which is the same formula. It agreed with
  `diversity.py` to its ten printed digits.
- `check_dadi.py` gives the folded projected spectrum with
  `dadi.Spectrum.from_data_dict(..., projections=[g], polarized=False)`.
  All 33 entries of the three populations agreed to ten digits. Run it
  with the 3.12 Python.
- `check_poppr.R` gives the alleles each population called, with
  `adegenet::genind2genpop`, the private ones, with
  `poppr::private_alleles`, and the expected heterozygosity per
  population, with `adegenet::Hs`. The first two agreed exactly. `Hs` is
  the plain expected heterozygosity with no correction for the sample,
  which is why it is not the reference for F_IS. It reads
  `panel_genotypes.tsv`, which `check_allel.py` writes.
- `check_allel.py` gives F_IS from `scikit-allel`'s per variant observed
  and expected heterozygosities, in the plain form, which is the form the
  spec checks against it, and writes `panel_genotypes.tsv`. Run it with
  `PYTHON_GIL=0` to keep the warning quiet.
- `fis_sim.py` is the measurement that chose between the two ways of
  summarizing F_IS over the variants. It draws genotypes at a known F,
  with the allele frequencies from the neutral spectrum, and reports the
  bias and the spread of the ratio of the two means and of the mean of the
  per variant ratios, with the plain and the unbiased expected
  heterozygosity, over 12 settings of the number of individuals, the
  rarest allele and the true F. The ratio of means with the unbiased
  expected heterozygosity was never more than 0.0056 from the truth and
  the mean of ratios was up to 0.0602 from it, which is the measurement
  the spec quotes.
