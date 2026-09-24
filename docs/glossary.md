# The names popnei uses

September 2026. popnei is a population genetics library in Rust, used from
Python and from TypeScript, the successor of the Python library pyNei. Its specs, its plans and
its code are written by different people and sessions, and this document
fixes which word each of them uses for the things of the domain, so that a
thing has one name in the docs, in Rust, in the binding and in Python. It
decides names and nothing else. How the variants flow is in
`architecture.md`, and what a calculation does is in its spec under
`docs/specs/`.

A document still says what a term means where it first uses it, as the
`writing` skill asks, because its reader may have only that page.

Each entry has the name, what it is, the name pyNei gives it when that
differs or when it is an identifier to look for, and the words that are not
used for it. The pyNei names were read in pyNei at commit ef0ca6e.

## The data

**individual.** One organism that was genotyped, known by its name. It is
the word of population genetics. VCF calls it a sample, and it is a
column of a VCF, and pyNei calls it a sample too, `samples`,
`num_samples`, `filter_samples`. In popnei the word is individual in the
prose and in the identifiers, `individuals`, `num_individuals`, and
"sample" is written only for a name of pyNei and for a statistical
sample. Inside the core an
individual is an index into the individuals of the reader. Not used:
sample, accession.

**population.** A named set of individuals that a calculation treats as a
group. "pop" in identifiers, `pop`, `pops`, and "population" in prose.
`pops` is a dict of population name to the names of its individuals in
Python, as in pyNei, and the indices of the individuals of each
population inside the core. Not
used: group, deme.

**variant.** One site of the genome with its alleles and the genotype of
every individual at it, a data line of a VCF. It can have more than two
alleles. Not used: SNP, which is one kind of variant, marker, locus.

**allele.** One of the forms a variant has, held as an `i8`: 0 is the
reference allele and 1 and above are the alternative ones, in the order of
the VCF.

**missing allele.** An allele that was not called, `.` in a VCF. It is the
constant `MISSING_ALLELE`, -1, in pyNei and in popnei, and never a literal
-1.

**ploidy.** How many alleles the genotype of one individual holds, the
same for every individual and every variant of a dataset.

**genotype.** The `ploidy` alleles of one individual at one variant, each of
which can be missing on its own. `gts` in identifiers.

**called genotype**, **missing genotype** and **half called genotype.** A
genotype is called when none of its alleles is missing, and missing when
at least one is, which is what `_calc_gt_is_missing` of pyNei's
`gt_counts.py` computes. A half called genotype, `0/.` in a VCF, has some
alleles called and some missing, so it is a missing genotype: it adds to
the missing rate of its variant, the missing genotypes divided by the
individuals, and it is not among the genotypes of the observed
heterozygosity. The calculations that count alleles, the
allele frequencies and the expected heterozygosity, still count the
alleles of it that were called, and their spec says so. Not used: partial
genotype, no call.

**observed heterozygosity.** Of a variant over some individuals, the
heterozygous genotypes divided by the called genotypes, where a genotype
is heterozygous when it is called and its alleles are not all the same.
`obs_het` in identifiers.

**called alleles.** How many alleles of a population at a variant are not
missing, the denominator of its allele frequencies.

**expected heterozygosity.** Of a variant in a population, one minus the
sum over the alleles of the frequency of each to the power of the ploidy:
the chance that as many gene copies as a genotype holds, taken at random
from the population, are not all alike. The unbiased one corrects it for
the frequencies being estimated from the copies it is computed over, Nei's
c/(c - 1) at ploidy 2 with c the called alleles, and is a statistic of its
own.
`exp_het` and `unbiased_exp_het` in identifiers. `docs/specs/stats.md`.

**polymorphic variant.** In a population, a variant whose major allele
frequency is below the polymorphism threshold, 0.95 by default, strictly;
a variable one has it below 1. `poly` in identifiers, as in pyNei's
`poly_threshold` and `num_poly`.

**major allele.** The allele of a variant with the highest frequency among
the called alleles of the individuals considered. How a tie is broken is for
the spec of the calculation to say.

**maf.** In pyNei and in popnei, the frequency of the major allele, as in
`filter_by_maf` and its `max_allowed_maf`. In most of the literature and in
plink2 the same three letters are the frequency of the minor allele, so a
text writes "the major allele frequency" in full where it first uses it.

**dosage.** For one genotype, how many of its alleles are not the major
allele of the variant, from 0 to the ploidy, every allele other than the
major one counting the same. The dosage matrix is the variants x
individuals array of them. pyNei: `to_012` and "the 012 matrix".

**fall-off curve.** The r² that two variants of a population are expected
to be in, against the recombination between them: Hill and Weir (1988)
with the correction of Weir and Hill (1986) for r² being measured on a
sample of individuals and not on the whole population. popnei fits it to
the pairs of each population, and the one number it fits is the ρ per
base pair below, the individuals of the population entering it as they
are. `decay` in identifiers, as in `fit_ld_decay` and `decay_per_pop`,
and "the fall-off curve" or "the fitted curve" in prose.
`docs/specs/ld.md`. Not used: the LD decay curve, the decay model.

**ρ per base pair.** The one number popnei fits to the fall-off of r²
with distance: 4Nr, four times the effective size of the population times
the recombination per base pair, so that two variants d base pairs apart
are separated by a scaled recombination ρ of d times it. `rho_per_bp` in
identifiers, and "the ρ per base pair" in prose. popnei cannot tell the
effective size and the recombination apart and gives their product.
`docs/specs/ld.md`. Not used: C, which is the letter the papers give it,
and decay rate.

**r² at distance 0.** The value of the fall-off curve where the two
variants are 0 base pairs apart, its own ceiling: two variants that never
recombine still do not reach an r² of 1, because their allele frequencies
drift apart. How many individuals the population has fixes it on its own,
0.46198347107438015 at 100 of them, and no pair of the dataset moves it.
It is what the half distance below is half of. `r2_at_zero` in
identifiers. `docs/specs/ld.md`. Not used: intercept, plateau, which is
what a reader may call the 1 over the individuals the curve falls towards
instead.

**half distance.** The distance at which the fall-off curve above has
fallen to half of its r² at distance 0. It is the fall-off of linkage
disequilibrium of a population as one number, and the one a web
application plots. It is read off the curve, so it can be beyond every
distance the fit was given. `half_dist` in identifiers.
`docs/specs/ld.md`. Not used: half life, LD decay distance.

**component.** A principal component: one of the directions, at right
angles to each other, along which the individuals of a standardized table
vary most, the first the one with the largest variance. "PC" in the names
of a result, `PC0`, and `comps` in identifiers. Not used: axis, eigenvector,
which is how a component is computed and not what it is.

**projection.** Where an individual falls along a component, the
coordinate a user plots. pyNei: `projections`. Not used: score, which is
R's word, coordinate.

**princomps.** The weights of each variant, or of each trait, in each
component, components x variants, a field of the result of a PCA under
the name pyNei gives it. "Weight" in prose. Not used: loading, rotation,
which is R's word.

**Kosman distance.** The distance between two individuals of Kosman and
Leonard (2005): at a variant, the alleles of the two genotypes that do
not pair with an equal allele of the other, over the ploidy, which for
diploids is 0 for the same genotype, 1 for two genotypes with no allele
in common and 0.5 otherwise, averaged over the variants at which both
genotypes are called. `docs/specs/dists.md`.

**Hudson's F_ST.** The distance between two populations that says how
much of the diversity the two hold together lies between them rather than
within them, estimated as Bhatia et al. (2013) recommend for SNPs: the
sum over the variants of the between population heterozygosity minus the
within one, over the sum of the between one. `fst` in identifiers.
`docs/specs/dists.md`. Not used: the fixation index, and F_ST alone where
a text could mean Nei's G_ST, which popnei also gives and which is a
different number.

**f_2.** The distance between two populations that is how far their allele
frequencies have drifted apart, the numerator of Hudson's F_ST over the
variants that counted, with the sampling bias taken out. It adds up along
a tree, which is what f_3 and f_4 are built on. `f2` in identifiers.
`docs/specs/dists.md`.

**chord distance.** The distance between two populations of Cavalli-Sforza
and Edwards: the square root of every allele frequency puts each
population on a sphere of radius 1, and the distance is the straight line
between them. It is Euclidean, so a principal coordinate analysis of a
matrix of them has no negative eigenvalues. Its square is Nei's D_A.
`chord` and `da` in identifiers. `docs/specs/dists.md`.

**Jost's D.** The distance between two populations that says how much of
their allelic variety is not shared, 0 when they have the same alleles at
the same frequencies and 1 when they share none. It answers a different
question from F_ST, and it is the one to read on markers with many
alleles. pyNei: `calc_jost_dest_pop_dists`, and `dest` in identifiers,
after the D_est the literature writes. `docs/specs/dists.md`.

**G_ST.** Nei's fixation measure between two populations, the share of the
diversity of the two that lies between them, from the expected
heterozygosities corrected for the sample as Nei and Chesser do. It cannot
reach 1 when the populations are diverse: with two of them its ceiling is
(1 - H_S)/(1 + H_S). `gst` in identifiers. `docs/specs/dists.md`. Not used:
F_ST for it, which in popnei is Hudson's and a different number.

**G''_ST.** G_ST rescaled so that it reaches 1 when the two populations
share no allele, whatever their diversity, as Meirmans and Hedrick (2011)
define it. `gst_standardized` in identifiers, since the literature's name
is not an identifier. Hedrick's earlier G'_ST, which divides G_ST by its
ceiling, is a different number and popnei does not give it.
`docs/specs/dists.md`.

**resampling group.** The variants that a standard error leaves out
together: a stretch of one chromosome, or one variant. The literature
calls it a block and calls the method the block jackknife; popnei says
group, because a block here is the run of variants a reader gives.
`jackknife_group` in identifiers. `docs/specs/dists.md`. Not used: block,
window.

**distance vector.** The distances of every pair of N individuals or
populations as one array, in the order (0, 1), (0, 2), ..., (0, N-1),
(1, 2), ..., the upper triangle of the square matrix row by row.
`dist_vector` in identifiers, as in pyNei. Not used: condensed matrix,
the name scipy gives the same order.

**standardized dosage.** The dosage of a genotype with the mean dosage of
its variant taken from it and the result divided by
`sqrt(ploidy * p * (1 - p))`, where p is the mean dosage over the ploidy:
the spread the allele frequency of the variant gives it under Hardy
Weinberg. It is what the kinship is built from, and it is not the dosage
divided by its own standard deviation, which is what a PCA of the variants
standardizes with. `z` in the formulas, as in the literature.

**kinship.** The genomic relationship matrix of VanRaden (2008), which GCTA
and plink2's `--make-rel` also compute: for every pair of individuals, the
standardized dosages of the two multiplied together and summed over the
variants, divided by the per pair denominator. An entry off the diagonal is
twice the coancestry of the pair and one on the diagonal is 1 plus the
inbreeding of that individual. `docs/specs/kinship.md`. Not used: GRM,
relationship matrix, K, which is what the formulas call it.

**per pair denominator.** How many variants have a called genotype in both
individuals of a pair, which is what that pair's entry of the kinship is
divided by. With no missing genotype it is the same number for every pair.
`num_vars_per_pair` in pyNei. `docs/specs/kinship.md`.

**trait.** What a user measured on each individual and wants the variants
tested against: **continuous**, a measurement, or **binomial**, 0 or 1.
`trait` in the arguments and `TraitType` in the types. The value itself,
one number per individual, is the **phenotype**, as in pyNei, which is the
argument a user passes. `docs/specs/gwas.md`.

**covariate.** A number per individual whose effect on the trait has to be
taken out but is not what is being tested, such as the sex or the field a
plant grew in. `covariates` in the arguments, a frame indexed by
individual. Not used: fixed effect, which is R's word and which in a mixed
model also covers the variant.

**design.** The matrix of one row per tested individual and one column per
number a model fits: a column of ones for the intercept and one for each
covariate. `design` in the core crate, and `d` in the formulas.

**null model.** The model of the trait fitted once with the covariates and
the kinship in it and no variant, which every variant is then tested
against. `NullModel` in the results. `docs/specs/gwas.md`.

**mixed model.** A model with the kinship in it as the covariance of a
random effect, so that related individuals are expected to resemble each
other before any variant is looked at. The two of popnei are the linear
mixed model, `lmm`, and the logistic one, `glmm`, two of the four values of
`GWASModel`.

**Wald test.** The test of a variant that fits the model again with the
variant in it and asks how many of its own standard errors the effect is
away from 0. `TestType.WALD`.

**score test.** The test of a variant that never fits the model with the
variant in it: it asks how steeply the fit would improve if the effect were
let off 0, measured at the null model. `TestType.SCORE`.

**heritability.** The variance of the kinship effect over the sum of it and
the residual variance: the share of the trait's variance the kinship
explains. A field of `NullModel`, and only the linear mixed model has one.

## How the data moves

**block.** Consecutive variants held as contiguous arrays, the `Block`
struct, which is how the variants flow from a source to a calculation: a
reader gives blocks, a filter compacts them and a calculation consumes
them, as rows or as a matrix. One variant of a block is a view into it,
`VariantRef`, and "row" is that variant as a line of the arrays. A block holds
about 5 million genotypes, the size pyNei gives its chunks, which is a few
thousand variants.
pyNei: chunk, `VariantsChunk`. "Chunk" is used only for pyNei's own. What
the Python `Variants` of popnei gives from `iter_blocks` is a block. Not
used: batch, which is written for two other things, arrow's unit of the
vars file, and the lines that the VCF reader reads and parses together,
several to a block; and window, which has a meaning of its own below. A
vars file is written with one batch for each block, and read back in
blocks of any size.

**window.** The variants that a calculation or a filter compares one
variant with: those on its chromosome whose position is no more than a
stated distance from it. It is never a run of blocks nor a number of
variants: a window is a stretch of a chromosome in base pairs, and how
many variants fall in it is whatever the dataset has there. The filter by
linkage disequilibrium of `docs/specs/filters.md` holds the variants it
has kept inside the window of the variant it is looking at, and the curve
of linkage disequilibrium against distance of `docs/specs/ld.md` compares
each variant with the ones inside its own. Not used: window for a block
or for a run of blocks, which is what a reader holds and not what a
calculation compares.

**member.** One gzip stream of a gzipped file. A file that bgzip wrote is
many of them one after another, each with 64 KiB of text at most and each
stating its own size, and the reader of a VCF cuts them by those sizes and
decompresses one at a time, as `docs/specs/io_vcf.md` describes.
BGZF: block. "Block" in popnei is a block of variants, so the word for a
gzip stream of such a file is member, which is what the gzip format calls
it, and the empty member at the end of a bgzipped file is the mark of its
end and not the empty block.

**region.** A stretch of one chromosome, from a smallest to a largest
position, both included. The vars file keeps, for each of its batches, the
region of every chromosome that has variants in it.

**reader.** Anything that gives blocks through the `BlockReader` trait: the
VCF reader, the vars file reader, and a filter or `reblock`, which are
readers over another reader. Not used: record, record level and block
level, the words of the first version of `architecture.md`, which had a
single variant that a reader filled.

**vars file.** The arrow file in which popnei keeps variants, feather v2,
which any program with an arrow library opens as a table. Its format is
popnei's own, in `docs/specs/io_vars.md`. pyNei has a file of the same
name and another format, and neither library reads the other's.

## The layers

popnei is used from two languages, Python and TypeScript. Both reach the
same core crate, each through a binding crate with a package on top, as
section 8 of `architecture.md` lays out.

**core crate.** `crates/popnei`, pure Rust with no Python and no
JavaScript in it, where every calculation is.

**binding crate.** A crate that translates between another language and
the core crate and holds no calculation. There are two, and a text that
**step.** One entry of the list that a `Variants` holds besides its
source, a filter with its threshold or the filter of individuals with
its names. A step is added with a method of the
`Variants`, which returns nothing, it is run in every pass that starts
after it was added, and `variants.steps` lists them.

**consumer.** What runs a `Variants`: the function of a calculation,
`write_vars`, the method `iter_blocks`. It does not change the `Variants`,
and what it returns has the pass stats.

**pass stats** and **filtering stats.** What a pass counted, `PassStats`,
which every result of a consumer has as `pass_stats`: how many variants
the consumer took, and the filtering stats, `FilteringStats`, how many
variants each filter was given and how many it kept. pyNei keeps the
filtering stats in its `Variants` and gives them with
`gather_filtering_stats`.

**pass.** One reading of a source of variants from its start to its end,
through the steps that its `Variants` had when it started. A consumer
makes as many as its algorithm needs, each with readers and filters of its
own.

could mean either says which. "The binding crate" alone is used for what
holds for both.

**Python binding crate.** `crates/popnei-python`, written with pyo3. Its
Python module is `popnei._core`.

**linalg crate.** `crates/popnei-linalg`, the linear algebra of popnei,
its products and decompositions, with two backends behind one interface:
BLAS and LAPACK natively, and faer, a library written in Rust, in wasm,
where there is no BLAS, and natively too when the cargo feature `blas`
of the crate is off. It is the one crate of popnei with `unsafe` in it,
the calls to BLAS and LAPACK. The core crate calls it.
`docs/specs/linalg.md`. Not used: backend for the crate itself, which is
the word for each of its two libraries.

**JavaScript binding crate.** `crates/popnei-js`, written with
wasm-bindgen, the Rust tool that generates the JavaScript that calls the
exported Rust functions and their TypeScript declarations.

**Python package.** `python/popnei`, what a Python user imports: the
functions with pyNei's signatures, the result objects and the Python
`Variants`.

**TypeScript package.** `js/popnei`, what a web application imports: the
functions of the Python API with their names in camelCase, and the result
objects with typed arrays.

**pyodide wheel.** The Python binding crate and the Python package built
for pyodide, the Python that runs in a browser tab, for the person who
writes Python in a notebook there. Compiled for
`wasm32-unknown-emscripten`.

**wasm package.** The npm package that a web application installs: the
TypeScript package with the JavaScript binding crate compiled for
`wasm32-unknown-unknown` inside. No Python runs in the tab. Not used:
"the wasm build" when the text could mean either this or the pyodide
wheel.

## How results are verified

**reference program.** A program outside the project whose output a
calculation of popnei is checked against: plink2, GMMAT, rrBLUP, R. Not
used: tool, gold standard.

**literal.** A number from a reference program or from pyNei that is
written into a test as it was printed, with the spec saying how it was
got. A test compares against literals and does not compute what it
expects.

## When a name is missing or wrong

Where pyNei has a name for a thing that a user sees, a function, an
argument, a field of a result, that is the name. A thing of the domain
that several documents will speak of and that has no entry here gets one
in the commit that first needs it. An entry that turns out wrong is
changed here first, and then wherever a search for the old word finds it.
