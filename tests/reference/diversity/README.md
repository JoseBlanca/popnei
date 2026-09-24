# The reference numbers of the diversity module

24 September 2026. This directory holds the numbers that the tests of the
`diversity` module of popnei compare against, and the three scripts that
write them. The module gives, for each population of a dataset, the alleles
it called, the private alleles among them, the variants that vary in it,
those three also taken down to a common number of called alleles, the
folded site frequency spectrum at that number, and F_IS.
`docs/specs/diversity.md` says what each of the five is, and its "How it is
verified" of each item says which of the six files below checks it. pyNei
has none of the five, so nothing here is a comparison with pyNei: five of
the six files hold what a program outside popnei gave, and the sixth holds
values computed twice in exact arithmetic, because no outside program
computes them at all.

Two terms are used throughout and both are arguments of the module.
`min_num_individuals`, 20 everywhere here, is how many called genotypes a
population needs at a variant for the variant to count for it.
`num_called_alleles`, 20 in the two scripts that read the panel and
smaller in the made up cases of the third, is the number of called alleles
every population is brought down to: draw that many of the called alleles
a population has at a variant, without replacement, and take the
expectation over every such draw, which makes a population of 48
individuals and one of 84 comparable. A value taken that way is called
standardized below, and a draw of 20 is that operation.

## The panel the first two scripts read

`make_reference.R` and `make_reference.py` read the panel of
`docs/specs/stats.md`, `tests/reference/stats/panel.vcf.gz` with the
population of each individual in
`tests/reference/stats/panel_pops_bcftools.txt`: 1200 biallelic diploid
variants of 200 individuals, 3 in 100 genotypes missing whole, in the
populations `p0`, `p1` and `p2` of 48, 68 and 84 individuals. Both take the
allele counts of each population at each variant out of that VCF
themselves, and neither reads a file the other writes. At
`min_num_individuals` 20 all 1200 variants count for all three
populations, and every population has at least 20 called alleles at every
variant, so the outside programs, none of which has such a threshold,
count the variants popnei counts. Each script stops if a population of the
panel falls below either of those two numbers at any variant.

`enumerate_private.py` reads no file at all.

Each script compares what it got with the literals of the tables of
`docs/specs/diversity.md` and stops at the first one that differs, so the
files here hold the numbers of the spec or the script fails.

## make_reference.R: adegenet, poppr and vegan

    Rscript tests/reference/diversity/make_reference.R

It runs under R 4.6.1, at `/opt/homebrew/bin/Rscript` on this machine, and
stops on any version of its three packages other than `vegan` 2.7.6,
`adegenet` 2.1.11 and `poppr` 2.9.8. It takes 1.4 s on the machine named
at the end of this page. It writes three files, one per item of the spec,
and the name of a column says which program gave the numbers in it.

`panel_num_alleles.tsv`, one row per population:

- `num_vars_with_data`, the variants that counted for it, 1200 for each.
- `total_adegenet`, the alleles it called over those variants, from
  `adegenet::genind2genpop`: 2373, 2377 and 2384.
- `mean`, that total over the variants.
- `num_vars_in_draw`, the variants at which it has at least 20 called
  alleles, 1200 for each.
- `in_draw_vegan`, the alleles a draw of 20 is expected to show, averaged
  over those variants, from `vegan::rarefy(counts, sample = 20)`:
  1.9283948650, 1.9219209943 and 1.9197370844.

`panel_private_alleles.tsv`, one row per population: `num_vars_every_pop`,
the 1200 variants every population has data at, `total_poppr`, the alleles
the population called that no other population called, from
`poppr::private_alleles(count.alleles = FALSE)` called on the `adegenet`
genotype object that the script builds out of the VCF, which is 0, 0 and
1, and `mean`, that total over the variants.

`panel_variable_vars.tsv`, one row per population: `total_adegenet`, the
variants that vary in it, 1173, 1177 and 1184 of 1200, `ratio`, that total
over the variants, and `in_draw_vegan_minus_one`, the chance that a draw
of 20 varies, averaged over the variants. That column is named as it is
because it is not a second measurement: it is the `in_draw_vegan` of
`panel_num_alleles.tsv` minus 1. On a variant of two alleles a draw shows
one allele or two, so the alleles it is expected to show are one plus the
chance that it varies. `vegan` measures one of the two numbers and the
other follows. The script stops if any variant of the panel has more than
two alleles, which is what that identity needs, and on every run it checks
the identity itself on the 14 pairs of allele counts and draw size that
"How it is verified" of the variable variants of the spec gives.

Every float in the three files is written to 17 significant digits and not
to the ten the spec prints, because the tests compare with `vegan` within
1e-12 of the value and a number rounded to ten decimals is up to 5e-11
away from it.

## make_reference.py: dadi and scikit-allel, each in an environment of its own

    uv run python tests/reference/diversity/make_reference.py

It writes two files.

`panel_folded_sfs_dadi.tsv` has one row per count of the rarer allele in a
draw of 20, `rarer_allele` 0 to 10, and one column per population. A value
is the variants expected to show that many copies of the rarer allele, and
each column sums to the 1200 variants that counted. It comes from `dadi`
2.4.4, `Spectrum.from_data_dict(dd, [pop], projections=[20],
polarized=False)`, read out of `fs.data`: `dadi` masks bin 0 and the bins
above 10 of a folded spectrum and popnei reports bin 0, which is a
difference of presentation and not of value.

`panel_fis_plain_allel.tsv` has one row per population and holds
`fis_plain`, from `scikit-allel` 1.3.13: one minus the ratio of the mean
observed heterozygosity to the mean expected heterozygosity, the
heterozygosities being the per variant `heterozygosity_observed` and
`heterozygosity_expected` of that package. The values are -0.0237536998,
-0.0258924700 and -0.0247472838. This is the plain form, with no
correction for the size of the sample. popnei returns the unbiased form,
which corrects each variant's expected heterozygosity by c / (c - 1) with
c the called alleles, so the stored value is not the value popnei gives:
"How it is verified" of "The inbreeding coefficient F_IS" of the spec gives
both and says how each is compared.

Both files carry every digit of the float, so that a test can compare
within a tolerance the file does not limit.

The two programs run in two Python 3.12 environments that the script makes
itself with `uv venv --no-project --python 3.12`, under
`tmp/diversity-reference/` of the repository, which `.gitignore` leaves
out. One holds `dadi` 2.4.4 and the other `scikit-allel` 1.3.13, and a
later run reuses an environment whose version is the one asked for and
makes it again otherwise. The same directory holds the program each of the
two processes runs, the allele counts `dadi` projects and the genotypes
`scikit-allel` reads, all written again on every run. With both
environments already there the script takes 0.7 s on the machine named at
the end of this page.

Neither program is a development dependency of `pyproject.toml`, and
nothing of popnei imports either: this script is their only caller and it
runs each in a Python of its own. Which build of Python 3.14 the project's
environment has is what decides whether they could be development
dependencies at all, and both builds are on this machine. On 3.14.5, the
build that keeps the global interpreter lock and the one `.python-version`
pins by its patch version, both install and import. On 3.14.7, the free
threading build, which `uv venv --python 3.14` picks by itself here and
which section 3 of `docs/rust_core.md` records as a trap of the wasm
build, `dadi` does not install: its `nlopt` dependency resolves there to
2.6.2, which is built from source and stops at `CMake must be installed to
build the following extensions: nlopt._nlopt`, where on 3.14.5 `nlopt`
resolves to 2.11.0 and comes as a wheel. On that same build, importing
`scikit-allel` turns the lock back on for the whole process, with a
warning that its compiled module `allel.opt.model` has not been marked as
safe to run without the lock. Both builds were tried on 24 September 2026.
The two environments of Python 3.12 keep the numbers here out of that
question.

## enumerate_private.py: no outside program at all

    uv run python tests/reference/diversity/enumerate_private.py

The standardized private alleles of a population are, at one variant, the
alleles expected to be in its own draw and in no other population's draw
of the same size. They are the one number of the spec that no program
outside popnei computes. popnei computes them with a closed form, the
formula of Kalinowski (2004), which for each allele multiplies the chance
that the population's own draw shows it by the chance that no other
population's draw shows it. ADZE, the program of Szpiech, Jakobsson and
Rosenberg (2008) and the only one that gives these values, evaluates that
same formula from that same paper, so agreeing with ADZE would say that
popnei copied the formula correctly and nothing about whether the formula
computes the sentence above.

So this script computes each value twice and compares the two. Once by
that closed form. Once by listing every draw each population can make,
weighting each draw by the number of ways it can be taken, running over
every combination of one draw per population, counting the alleles that
are in the population's own draw and in no other, and averaging. Both run
in `fractions.Fraction`, so they agree exactly or not at all.

What that comparison is worth, and what it is not. The enumeration never
writes the closed form down, so it catches any error in the algebra: the
decomposition into one term per allele, the chance that an allele shows in
a draw, and the product over the other populations. Dropping the factor
that asks for the allele to be in no other population's draw makes 19 of
the 22 pairs that agree differ, measured on 24 September 2026. What it
does not check is the one assumption the closed form makes, that the draws
of two populations are independent: the enumeration runs over the
combinations of one draw per population as a product measure, which is
that same assumption. For those 22 pairs, whose populations share no
individual, it is sound rather than circular, because draws from disjoint
sets of gene copies are independent as a fact of the sampling and not as
an assumption. Where two populations do share an individual the assumption is
false and an enumeration over allele counts is as wrong as the formula.

`enumerate_private.tsv` holds 23 pairs of a case and a population, one
line each, over ten cases. A case is one variant and not a dataset: it
gives what every population called at that variant, and the value of a
pair is the per variant value that popnei averages over the variants of a
population. Eighteen of the pairs are the ones "How it is verified" of
"The private alleles" of the spec lists: the three variants of the worked
example that have a draw of 4, and five cases made up for the check, of
two and three populations at draws of 2 and 3. Two cases are there because
those 18 leave parts of an implementation untested, nine of their values
being 0 or 1 and 14 of their 18 population slots having exactly 4 called
alleles: a single population, which gets every allele it called, and three
populations of 3, 5 and 6 called alleles, no two alike. The tenth case is
the shared individual below.

The ten fields of a line are `case`, its name, which for the eighteen is
the name the spec gives them; `enumerated_over`, which of the two
enumerations gave the line;
`allele_counts`, the copies of each allele that each population called,
the populations separated by ` | `; `num_called_alleles`, the size of the
draw; `population`, `pop1` for the first population of `allele_counts` and
so on; `closed_form` and `enumerated`, the two values to 17 decimal
places, which the cargo tests of the module assert and each of which reads
back as the float64 nearest the exact value, checked by the script and not
assumed; `closed_form_exact` and `enumerated_exact`, those same two as
exact rationals, `14/15` and not `0.9333333333`, for a reader checking a
line by hand; and `difference`, `closed_form` minus `enumerated` as an
exact rational.

One line of the file has a difference that is not 0, and it is there on
purpose. It is the case named "the shared individual, the one pair the
closed form gets wrong": two populations that are both the one diploid
individual `0/1`, at a draw of one allele. The closed form gives 1/2. The
two draws are draws of the same genotype and can never give an allele to
one population and not the other, so the truth is 0, and the line holds a
difference of 1/2. That pair is enumerated over the labelled gene copies
of the individuals rather than over allele counts, so that a copy two
populations share is drawn by both or by neither, which is what lets it
show the error; `enumerated_over` reads `labelled gene copies` there and
`allele counts` on the other 22. The script refuses to write the file
unless those 22 agree exactly and this one disagrees by exactly 1/2, so a
difference of 1/2 on this line is the file being right and any other
difference anywhere is the file not being written at all.

## What the spec checks and this directory does not hold

"The private alleles" of the spec also gives the standardized private
alleles of `p0`, `p1` and `p2` on the panel at a draw of 20,
0.0112196177, 0.0099715392 and 0.0089014974. No file here holds them and
no script here can compute them, because no outside program gives any
standardized private allele value. They come from
`docs/reports/diversity-method/panel.py`, which computes the five
quantities in Python as the spec defines them, and they are literals of a
pytest test, which therefore checks that popnei's Rust agrees with that
Python and nothing more. What checks the algebra of the formula is the
enumeration above, over the 22 pairs whose populations share no
individual, and what shows where the formula stops being right is the
23rd.

## hierfstat, which would have replaced two of the programs

`hierfstat` computes in one package both of the quantities that `vegan`
and `scikit-allel` are the references for here: its `allelic.richness` is
the same rarefaction as `vegan::rarefy`, and it gives F_IS. It is not
installed on this machine and does not install.
`install.packages("hierfstat")` on 24 September 2026 fetched `hierfstat`
0.5-11 and stopped at its dependency `RcppParallel` 6.2.1, whose
configuration ends with `error: RcppParallel requires cmake (>= 3.5);
cmake was not found`, and then at `gaston`, a second package `hierfstat`
needs, which needs `RcppParallel` too. `docs/specs/dists.md` recorded the
same failure on 23 September 2026. What is missing is `cmake`, the same
program `dadi` wanted on the free threading Python above, so a machine
that has it may well install all three; this one does not have it.

## The machine

The six files were written on 24 September 2026 on the owner's Apple M5
Pro under macOS 27.0, with R 4.6.1 holding `vegan` 2.7.6, `adegenet`
2.1.11 and `poppr` 2.9.8, and with `dadi` 2.4.4 and `scikit-allel` 1.3.13
in the two environments of Python 3.12.13 that `make_reference.py` makes.
Running the three scripts again on that day left all six files unchanged:
`git status --short tests/reference/diversity` printed nothing.
