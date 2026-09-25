# How r² is measured, and where the numbers of the ld spec come from

22 September 2026. The programs behind the tables of `docs/specs/ld.md`
and of the item "The filter by linkage disequilibrium" of
`docs/specs/filters.md`: what a missing genotype does to the correlation
between two variants, what the correct rule costs against pyNei's, which
variants the filter keeps, and what plink2 gives for all of it. r² is the
square of the correlation, across the individuals, between the dosages of
two variants, where the dosage of a genotype is how many of its alleles
are not the major allele of its variant; the specs say what each number
means and what was concluded, and this page says what each file is and
how to run it. Nothing here is part of popnei's build or of its tests.
The plan that builds the module took `make_ld.py` and the worked example
out of this directory on 22 September 2026: they are
`tests/reference/ld/make_reference.py` and `tests/reference/ld/example.vcf`
now, and `tests/reference/ld/run_plink2.sh` runs the plink2 commands below
for the two of them and compares what they give with what is stored. The
scripts here read the files those two write, so the paths below name them
where they are.

The machine was the owner's Apple M5 Pro, 18 cores, macOS 27.0, with
plink2 v2.0.0-a.7.7 M1, bcftools 1.24, R 4.6.1, and the Python
environment of popnei, numpy 2.5.3 on Accelerate and pyNei at commit
ef0ca6e.

Every script takes the directory it works in from the environment
variable `LD_WORK`, and uses the current directory when it is not set.
Run them with popnei's Python, `.venv/bin/python`, from the root of the
repository.

    export LD_WORK=/tmp/ld_work
    mkdir -p $LD_WORK
    cp tests/reference/ld/example.vcf $LD_WORK/

## The datasets

`tests/reference/ld/make_reference.py` writes `$LD_WORK/ld.vcf`, the
dataset the spec calls `tests/reference/ld/ld.vcf.gz`: two chromosomes of
250 biallelic variants each, 1000 bp apart, and 100 diploid individuals,
whose haplotypes come from four founders recombined along the chromosome
at a rate of 2 in 100 between one variant and the next, with 3 in 100
genotypes then set to missing. Its seed is `numpy.random.default_rng(29)`,
so it is the same file on every machine. 68 of its 500 variants end up
with one dosage in every called genotype.

It also writes `$LD_WORK/many.pynei.dosages.tsv`, the 500 x 50 dosages
pyNei's `to_012` gives for `tests/reference/vcf/many.vcf`, which a cargo
test compares every one of against what popnei counts.
`tests/reference/ld/run_plink2.sh` makes those two files and the two
plink2 matrices again, six files in all, and compares each with the copy
stored in git.

`tests/reference/ld/example.vcf` is the worked example of both specs: 5
variants of 6 individuals, one genotype missing and one variant of one
dosage, whose r² can be worked out by hand.

The repository already has the two other files the measurements use,
`tests/reference/dists/panel.vcf.gz`, 200 individuals and 1200 biallelic
variants with 3 in 100 genotypes missing, and
`tests/reference/vcf/many.vcf`, 500 variants of 50 individuals with 54
variants of more than two alleles and 257 half called genotypes.

## The matrices plink2 gives

`--r2-unphased` writes r²; `square bin` writes the whole matrix as
float64, row after row, where the text plink2 writes by default has six
digits. The four the scripts read:

    python tests/reference/ld/make_reference.py

    plink2 --vcf $LD_WORK/ld.vcf --double-id --allow-extra-chr \
           --r2-unphased square bin --out $LD_WORK/ld_r2

    plink2 --vcf $LD_WORK/example.vcf --double-id --allow-extra-chr \
           --r2-unphased square bin --out $LD_WORK/ex

    plink2 --vcf tests/reference/dists/panel.vcf.gz --double-id \
           --allow-extra-chr --r2-unphased square bin --out $LD_WORK/panel_r2

    plink2 --vcf tests/reference/vcf/many.vcf --double-id --allow-extra-chr \
           --vcf-half-call m --r2-unphased square bin --out $LD_WORK/many

The file plink2 writes is named for what it holds, `.unphased.vcor2.bin`
for r². Putting `--r-unphased` in the place of `--r2-unphased` gives the
correlation itself, with the sign that says whether the major alleles of
the two variants go together, in `.unphased.vcor1.bin`; that is what the
sentence of `docs/specs/ld.md` about the sign of variant 856 of the panel
is read from, and no script here reads it.

## What each script measures

`missing_rules.py` is the table of "Missing genotypes" of
`docs/specs/ld.md`: over the 1.4 million pairs off the diagonal of the
panel it compares plink2's r² with what three rules for a missing
genotype give, the individual left out of the pair, the missing genotype
given the mean dosage of its variant, and pyNei's, which leaves it in as
a dosage of -1. It prints one line per rule, the median, the 99th
percentile and the largest of the difference from plink2, and the median
r² of the panel beneath them. pyNei's `_calc_rogers_huff_r2` returns r
and not r², so the script squares what it gives. It needs `panel_r2`.

`half_called.py` is the measurement of **Open 2** of that spec: on
`many.vcf` it compares plink2's r² with popnei's two possible rules for
picking the major allele of a variant that has half called genotypes,
counting their called allele and not counting it. It needs `many`.

`bins.py` is the three tables of the curve of r² against distance: it
picks, for each population, the variants that pass its major allele
frequency, runs plink2 on those individuals and those variants alone, and
bins the pairs as the spec defines the bins. It needs `ld.vcf` and runs
plink2 itself.

`decay.py` and `decay.R` are the table of the curve fitted to the
fall-off of r² with distance, in the same "How it is verified", and the
numbers that "The curve that is fitted" and **Open 4** of that spec
quote. `decay.py` picks the variants and the individuals of each
population as `bins.py` does, runs plink2 on them and writes, for each
population, one line per distance with the pairs at that distance and the
sum of their r²; `decay.R` fits Hill and Weir's expected r² to those
lines twice, with R's `optimize` and with R's `nls`, and prints what the
two give, what the curve is at the middle of each of the ten bins, what
Sved's curve gives instead, what fitting ten and fifty bin means gives in
place of the pairs, what the three choices of **Open 4** give, what the
ends of the searched range do, the ρ at which the curve is half of its
value at 0 for four sample sizes, and the case of the first cargo test,
whose r² is taken from the curve itself. Run `decay.py` first: it writes
the three files `decay.R` reads. `decay.py` needs `ld.vcf` and runs
plink2 itself; `decay.R` needs R, 4.6.1 here.

`decay_truth.py` is the evidence behind the two choices that "The curve
that is fitted" of `docs/specs/ld.md` makes, which curve and what n is.
It runs a Wright-Fisher population of 200 individuals for 1000
generations, with two chromosomes of 1500 sites 200 bp apart,
recombination 2.5e-07 per bp and mutation 1e-04 per site, so that the
population reaches the drift and recombination the curve is the
expectation of and its 4Nr per bp is known to be 2e-04. It then fits Hill
and Weir's curve with n the individuals, with n the individuals times the
ploidy and with no sample term, and Sved's, to the pairs of one
chromosome, and prints how far each lands from that 4Nr and what each
gives at the shortest distance; and it takes the pairs of variants on
different chromosomes, which assort independently and so carry no linkage
disequilibrium, at four sample sizes, to measure whether the sampling
term goes by the individuals or by the gametes. Three seeds, about 40 s,
needs no file and no other program.

`prune.py` is the table of the filter: it walks the variants with the
rule of `docs/specs/filters.md` on plink2's own r² matrix, checks the two
properties the spec says hold of the set it keeps, and runs plink2's
`--indep-pairwise` beside it. It needs `ld_r2` and `ex`.

`residual_ld.py` and `threshold_knob.py` are the table of what popnei's
rule and plink2's `--indep-pairwise` cost each other and the sweep of the
threshold beside it, both in "How it is verified" of the filter item of
`docs/specs/filters.md`: how many variants each rule keeps and how much
r² is still standing between the variants it kept. They need `ld_r2`.

`speed.py` and `tiles.py` are the two tables of "Speed" of
`docs/specs/ld.md`, the cost of the six products against one product on a
block of 5000 variants and 1000 individuals, and the cost of a tile pair
of 256 and of 512 variants. They make their own dosages and need no file.

`pynei_first_variant.py` is the 1 of 6 against 5 of 6 under "What pyNei
does that is odd" of the filter item: pyNei's filter keeps the first
variant of the first chunk whatever it is, so when that variant has one
dosage every correlation against it is NaN and nothing else is ever kept.
It needs no file.
