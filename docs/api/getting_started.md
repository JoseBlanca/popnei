# Getting started

This page takes a VCF through the analyses most studies start with: the
statistics of each population, the F_ST between them, a principal component
analysis and a principal coordinate analysis of the distances between the
individuals. It also shows how to install popnei in a project of your own.
Every block of code on it was run on the panel that popnei's tests use,
`tests/reference/stats/panel.vcf.gz`: 1200 variants of 200 diploid
individuals in three populations, p0, p1 and p2, of 48, 68 and 84
individuals. The paths are relative to the root of the repository.

## Installing popnei in an analysis project

popnei is not on PyPI yet. It is installed from a wheel, the file that
holds the compiled package, which you build in the repository. Building
needs the Rust toolchain:

```console
$ cd ~/devel/popnei
$ uv run maturin build --release
```

The wheel is written to `target/wheels/`, and maturin prints its name,
which says the Python and the platform it is built for. The name below is
the one for a Mac with an Apple chip. In the directory of the analysis,
make a project, pin Python 3.14.5 and add the wheel:

```console
$ cd ~/analyses/my_study
$ uv init --bare
$ uv python pin 3.14.5
$ uv add ~/devel/popnei/target/wheels/popnei-0.1.0-cp314-cp314-macosx_11_0_arm64.whl
```

The pin to 3.14.5 is needed. With a bare `3.14`, uv takes the newest 3.14
interpreter it has, which on a Mac can be the free-threaded build. The wheel
is not built for that interpreter, and `uv add` then fails with "the binary
distribution is incompatible with the current platform". Adding the
repository by its path instead of by the wheel does not work either: uv
installs numpy and pandas and leaves popnei out, because the repository's
`pyproject.toml` tells uv not to install popnei into its own development
environment.

After a change to popnei, build the wheel again and run
`uv add --reinstall-package popnei <wheel>` in the analysis project.

## Opening a VCF

{func}`popnei.open_vcf` reads the header of the file and gives back a
{class}`popnei.Variants`: a handle over the file, with the options it is
read with. It holds no genotypes. Each calculation you give it reads the
file again, from its start, and one reading of the file from start to end
is called a pass.

```python
import pandas
import popnei

variants = popnei.open_vcf("tests/reference/stats/panel.vcf.gz")
variants.num_individuals  # 200
variants.ploidy           # 2
```

A VCF may be plain or gzipped. `open_vcf` leaves out the variants whose
FILTER column is neither `PASS` nor a dot, unless it is given
`only_passed=False`. It reads a ploidy of 2 unless it is given another one,
and a genotype of any other ploidy is an error.

## The populations

A calculation that works per population takes the populations as a dict:
each key is the name of a population, and its value is the list of the
names of its individuals, as the header of the VCF writes them. The panel's
populations are in a table of two columns, the individual and its
population:

```python
table = pandas.read_csv("tests/reference/stats/panel_pops.txt", sep="\t")
pops = table.groupby("popcat")["IID"].apply(list).to_dict()
{pop: len(inds) for pop, inds in pops.items()}  # {'p0': 48, 'p1': 68, 'p2': 84}
```

## Filtering

A filter is a method of the `Variants`. It returns nothing and reads
nothing. It is added to the list of steps of the `Variants`, and every
later pass runs it inside the Rust core as it reads the file. The variant
filters are these five:

- {meth}`~popnei.Variants.filter_by_missing_data`, over the rate of missing
  genotypes of each variant
- {meth}`~popnei.Variants.filter_by_maf`, over the frequency of its major
  allele, the most frequent one; MAF here is the major allele frequency, as
  in pyNei, and not the minor one
- {meth}`~popnei.Variants.filter_by_obs_het`, over its observed
  heterozygosity
- {meth}`~popnei.Variants.filter_by_ld`, which drops the variants in
  linkage disequilibrium with one kept before them nearby
- {meth}`~popnei.Variants.filter_by_regions`, which keeps the variants in
  the regions of a BED file, or those outside them

{meth}`~popnei.Variants.filter_individuals` keeps the individuals you name
and drops the genotypes of the rest. Each of the first three filters keeps
the variants whose value is at most the threshold it is given. A `Variants` takes each of
the six filter methods once, and calling one of them a second time is an
error. To try a second threshold, call `open_vcf` again for a second
`Variants`, which reads the header only.

```python
variants.filter_by_missing_data(max_allowed_missing_rate=0.1)
variants.filter_by_maf(max_allowed_maf=0.95)
variants
# <Variants of tests/reference/stats/panel.vcf.gz, ploidy=2, only_passed=True,
#  missing_data(max_allowed_missing_rate=0.1), maf(max_allowed_maf=0.95)>
```

## Reading a result

Every calculation returns an object whose fields are pandas series and
frames, indexed by the names of the populations or of the individuals. Each
result also has a `pass_stats` field with the counts of its pass:
`num_vars` is how many variants the filters let through to the calculation,
and `filtering` how many each filter was given and how many it kept:

```python
diversity = popnei.calc_pop_diversity(variants, pops=pops)
diversity.pass_stats
# PassStats(num_vars=1175, filtering={
#   'missing_data': FilteringStats(vars_processed=1200, vars_kept=1200),
#   'maf': FilteringStats(vars_processed=1200, vars_kept=1175)})
```

Of the 1200 variants of the panel, none has a missing rate above 0.1, and 25
have a major allele frequency above 0.95.

## Statistics of each population

{func}`popnei.calc_pop_diversity` counts the alleles of each population,
its private alleles and its variable variants, and calculates its F_IS.
{func}`popnei.calc_per_var_distribs` calculates the statistics of each
variant (observed and expected heterozygosity, major allele frequency,
missing rate, and whether the variant is polymorphic) and gives, for each population, their mean and their
histogram.

```python
diversity.fis
# p0   -0.013206
# p1   -0.017993
# p2   -0.018406

distribs = popnei.calc_per_var_distribs(variants, pops=pops)
distribs.exp_het.mean
# p0    0.353492
# p1    0.353588
# p2    0.348928
```

A variant counts for a population only when at least 20 of its individuals
have a called genotype at it. The argument `min_num_individuals` changes
that number.

## Distances between populations

{func}`popnei.calc_pop_dists` gives seven measures for every pair of
populations: Hudson's F_ST (`fst`), f_2 (`f2`), the chord distance of
Cavalli-Sforza and Edwards (`chord`), Nei's D_A (`da`), Jost's D (`dest`),
Nei's G_ST (`gst`), and the G''_ST of Meirmans and Hedrick
(`gst_standardized`), G_ST rescaled to reach 1 for two populations that
share no allele.
The reference page describes each one.

The standard errors are calculated by jackknife: the variants are cut into
groups, each group is left out in turn, and the spread of the estimates
without it gives the error. The argument `jackknife_group` says how the
groups are cut, and it has no default. It is a length in base pairs of a
chromosome, so that variants in linkage with each other are left out
together, `"variant"` for one group per variant, or `None` for no standard
errors. Each measure is a {class}`popnei.Distances`, and
its `square_dists` is the square matrix as a frame:

```python
pop_dists = popnei.calc_pop_dists(variants, pops, jackknife_group=None)
pop_dists.fst.square_dists
#           p0        p1        p2
# p0  0.000000  0.105259  0.102994
# p1  0.105259  0.000000  0.109913
# p2  0.102994  0.109913  0.000000
```

## Principal component analysis

{func}`popnei.do_pca_from_variants` makes each variant one column of the
table it analyses. The value of an individual in that column is its dosage:
how many alleles of its genotype are not the major allele of the variant,
0, 1 or 2 for a diploid. A variant with more than two alleles is an error
unless the function is given `transform_to_biallelic=True`, which counts
every allele other than the major one alike. The function gives the first
10 components unless it is given `num_prin_comps`:

```python
pca = popnei.do_pca_from_variants(variants)
pca.explained_variance_percent.head(3)
# PC000    7.725980
# PC001    5.607519
# PC002    1.562752
pca.projections.iloc[:3, :3]
#          PC000      PC001     PC002
# s000  1.533907  12.930970 -5.957159
# s001  2.802061  13.383321 -8.130564
# s002  2.307825  12.647294 -6.461064
```

`projections` has one row per individual and one column per component,
ready to plot.

## Distances between individuals and principal coordinates

{func}`popnei.calc_pairwise_kosman_dists` gives the Kosman distance of every
pair of individuals. At one variant it is the share of the alleles of the
two genotypes that differ, 0, 0.5 or 1 for diploids, and the distance of the
pair is its mean over the variants at which both are called. A principal coordinate analysis places the individuals on axes
from these distances, as a principal component analysis does from a table.

The distances of the panel are not Euclidean: no set of points in any space
has them as its distances. {func}`popnei.do_pcoa` refuses such distances
and says how far from Euclidean they are:

```python
dists = popnei.calc_pairwise_kosman_dists(variants)
popnei.do_pcoa(dists)
# ValueError: 45 of the 200 eigenvalues of the matrix of the squared
# distances are negative, 3.02 percent of the sum of all of them, so the
# distances are not Euclidean ...
```

Lingoes' correction adds the same constant to the square of every distance,
the smallest one that makes them Euclidean.
{func}`popnei.correct_dists_by_lingoes` applies it:

```python
corrected = popnei.correct_dists_by_lingoes(dists)
pcoa = popnei.do_pcoa(corrected.dists)
pcoa.explained_variance_percent.head(2)
# PC000    9.651062
# PC001    6.643066
```

{func}`popnei.do_pcoa_from_variants` does the distances, the correction and the analysis in one pass over the
variants when it is given `correct_by_lingoes=True`.

## What is not on this page

The reference pages have the rest: the vars file, popnei's own format for a
dataset that is read many times ({func}`popnei.write_vars` and
{func}`popnei.open_vars`); writing the filtered variants back to a VCF
({func}`popnei.write_vcf`); the kinship; the decay of linkage
disequilibrium with distance; and the association studies. Each docstring
says how the function differs from its pyNei counterpart.
