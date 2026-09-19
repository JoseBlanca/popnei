# The stats module: expected heterozygosity per population

September 2026. The expected heterozygosity tells a user of popnei how much
genetic diversity each of their populations has, in the plain form and in the
form corrected for the number of samples. There is no code. This part of the
spec of the `stats` module covers that one calculation, the `expected het` of
the module's row in section 9 of `docs/architecture.md`. It uses two things
written elsewhere: the row helper of the `variant` module that counts how often
each allele appears in one variant over a set of sample indices, and the item
of this spec for `calc_per_var_distribs`, the one pass that calculates the per
variant statistics and summarizes each of them into a mean and a histogram.
That item and the others of the module, the allele counts and frequencies per
population, the observed heterozygosity, the per sample statistics and the
polymorphism ratio, are not written yet.

## Expected heterozygosity

### What it gives

For one variant and one population, the chance that gene copies taken at
random from that population are not all of the same allele: the genetic
diversity of that population at that site, 0 when every copy carries the same
allele and close to 1 when many alleles are at similar frequencies. Over a
dataset the user gets its mean over the variants and its histogram.

The plain one is

    H = 1 - sum over a of p_a^k

where a runs over the alleles of the variant, p_a is the frequency of allele a
among the alleles of the population that were called at that variant, and k is
the ploidy. For k = 2 this is the classic 1 - p² - q², the chance that two
copies differ. For any other ploidy it is the chance that k copies are not all
alike, and not the chance that two of them differ.

Two more counts are needed for the unbiased one. c is how many alleles of the
population were called at that variant, the denominator of every p_a. n is
c / k, how much called data there is counted in genotypes, and it can be a half
when a genotype is half called. The unbiased one is

    H_u = (2n / (2n - 1)) * H

At k = 2, where 2n = c, that factor is c / (c - 1). The frequencies p_a are
estimated from the same c copies the statistic is computed over, and that makes
H too small on average, in the way that dividing a sum of squares by c rather
than by c - 1 makes a variance too small; multiplying by c / (c - 1) takes the
bias out. It is Nei's 1978 correction. GenAlEx, the population genetics add-in
for Excel, prints it as the unbiased heterozygosity for codominant data, the
name pyNei's docstring gives it. popnei keeps it as the default, as pyNei does,
so that a user comparing against GenAlEx or against pyNei sees the same number.

Every allele of the variant counts. A multiallelic variant is not collapsed to
the major allele against the rest here, as it is in other calculations of
pyNei.

### What it is in Python

The statistic is one of the four that `calc_per_var_distribs` calculates in one
pass, asked for as `stats="exp_het"`, and what a user gets back is the field
`exp_het` of the result, a `StatsDistrib` with the mean of the variants that
have a value and their histogram, one column per population. The function, the
`pops`, `min_num_samples`, `hist_kwargs` and `num_threads` arguments it shares
with the other three, and `StatsDistrib`, are the item of this spec for that
function. Two arguments are this statistic's own: `unbiased_exp_het=True`,
which chooses between the two formulas, and `ploidy=None`, which takes the
ploidy of the variants when it is not given. It mirrors `_calc_exp_het_per_var`
and `_calc_unbiased_exp_het_per_var` in `pynei/diversity.py`, which get their
frequencies from `_count_alleles_per_var` in `pynei/gt_counts.py`.

The objectives ask for every difference from pyNei to be written down. One is
decided here. `_calc_exp_het_per_var` sorts the population names,
`sorted(pops.keys())`, while the observed heterozygosity and the major allele
frequency keep the order of the `pops` dict, so in pyNei the columns of the
three do not line up. popnei gives its results keyed by population name and the
Python layer puts them in the order of the `pops` dict for every statistic. No
value changes. The other differences are the open points.

### Missing genotypes, populations with little data, and what pyNei asserts

The frequencies are over alleles, not over genotypes, so a half called
genotype contributes its called allele. A population of five diploid samples
with four genotypes called and one half called has c = 9 and n = 4.5, while
the observed heterozygosity of the same population at the same variant counts
4 called genotypes. The two statistics disagree about how much data the
population has.

A variant has no value in a population when n < `min_num_samples`, strictly
less, so a population with exactly `min_num_samples` genotypes' worth of called
data keeps its value. The default is 20, inherited from pyNei's
`MIN_NUM_SAMPLES_FOR_POP_STAT`. Nobody has measured whether 20 is the right
threshold. With no value the variant is out of the numerator and the
denominator of the mean and falls in no histogram bin, so the mean is over the
variants that had enough data. A test on a small dataset has to lower
`min_num_samples` or every value is missing; pyNei's own tests pass 1 or 5.

A population with no called allele at a variant has no value, except in a chunk
where no sample is called at any variant (**Open 1**, below).

A diploid population with one called allele, which needs `min_num_samples` at
0 to get this far, has n = 0.5 and 2n - 1 = 0: pyNei gives 0 for the plain one
and NaN for the unbiased one. popnei gives `Some(0.0)` and `None`.

`test_calc_exp_het` in `test/test_diversity.py` asserts the per variant numbers
of the worked example below, over the same two populations, and the mean and
the histogram of the default population over all five samples.
`test_per_var_distribs_with_threads` in `test/test_threads.py` asserts that the
mean and the histogram counts do not change with 2 or 4 threads. No test of
pyNei varies the size of the chunks for this statistic, and the chunk with
nothing called below is a case where its mean would change with them. popnei
adds a test that the size of the block does not change the result, which holds
once Open 1 is answered as recommended.

### What pyNei does that is odd

The three runs below are of pyNei at commit ef0ca6e, each one call of the named
function from a python session on a chunk of `Variants.from_gt_array` written
for it.

`_count_alleles_per_var` takes the alleles it counts from the largest allele in
the whole chunk, `gts.max()`, because asking numpy which alleles are there cost
more than the counting. A chunk is pyNei's unit of work, a few thousand
variants, the block of section 2 of the architecture. When no genotype in the
chunk is called at all, the largest allele is the missing one and there are no
alleles to count: the sum over them is 0 and every variant of that chunk gets a
plain expected heterozygosity of 1 - 0 = 1 and an unbiased one of -0.0.
Measured on two variants of three samples with every genotype `./.`: plain 1.0
and 1.0, unbiased -0.0 and -0.0, and through `calc_per_var_distribs` a mean of
0.0 with both variants counted in the first histogram bin. The same variant in
a chunk where any sample is called anywhere gives no value (**Open 1**, below).

The unbiased factor is 2n/(2n-1) whatever the ploidy, in
`_calc_unbiased_exp_het_per_var`, although at ploidy k the population holds
c = kn copies and not 2n. Measured with `min_num_samples=1` on three tetraploid
samples, genotypes 0/1/2/3, 0/0/1/1 and 0/1/2/3: the unbiased value is 1.1685,
above 1, so it is in no bin of the default histogram while the mean counts it
(**Open 2**, below).

The `ploidy` argument, when it differs from the ploidy of the data, is used as
the exponent k and as the number of alleles the samples are expected to hold,
`len(pop) * ploidy - missing_alleles`, while the `min_num_samples` test keeps
the ploidy of the data, which is the one `_count_alleles_per_var` uses. Measured
with `ploidy=4` on four diploid samples with allele counts 5 and 3: a plain
0.8276 and an unbiased 0.9459, where the 8 alleles the samples really hold
would give 0.8276 · 4/3 = 1.1035 (**Open 3**, below).

### How it runs

At the record level. Per variant, the counts of each allele over the sample
indices of a population come from the row helper of the `variant` module, and
the value follows from those counts alone. What is kept from one variant to the
next is the accumulator of the `calc_per_var_distribs` item: per population a
sum, a count of the variants that had a value, and the histogram bins. The
memory does not grow with the variants or with the samples. pyNei's chunk wide
set of alleles has no counterpart here, because the sum runs over the alleles
the row itself has.

### How it is verified

Against plink2 v2.0.0-a.7.7, the arm64 build of 18 September 2026. Its
`--hardy` report prints, per variant, `E(HET_A1)`, the frequency of
heterozygotes expected under Hardy Weinberg, which for a biallelic variant is
1 - p² - q², the plain statistic. `--loop-cats` runs the report once per
category of a column of a file, which is how the populations are told apart:

    plink2 --vcf sim_missing.vcf --pheno cats.txt --loop-cats popcat \
           --hardy --nonfounders --out h

on the panel that `test/gwas_reference/make_reference.py` simulates, 200
samples, 1200 biallelic diploid variants, three subpopulations, 3 in 100 of the
genotypes missing whole; `cats.txt` holds `IID` and a `popcat` column, p0, p1
or p2, built from the `pop` column of `phenotypes.csv`. `--hardy` reports the
founders only unless `--nonfounders` is given, and every sample of a VCF is a
founder, so here the flag changes nothing: the two reports are identical files.

Against `_calc_exp_het_per_var(chunk, pops, min_num_samples=20)` on
`sim_missing.vars`, over the 1200 variants and the 3 populations, the largest
absolute difference is 5.0e-7 in each population, with no missing value on
either side, so the tests compare within 1e-6 absolute, half a unit of the last
of the six digits plink2 prints. The numbers written into the first of those
tests as literals are those of variant `var0000`: in p0, 39 homozygotes for the
reference, 8 heterozygotes, 1 homozygote for the alternative and an expected
heterozygosity of 0.186632; in p1, 15, 34, 18 and 0.498998; in p2, 81, 2, 0
and 0.0238061.

No program outside the project prints the unbiased one. What is checked instead
is that the factor rebuilds it from plink2's own numbers: c is twice the
genotypes it counted, and for `var0000` in p1, c = 134 and
0.498998 · 134/133 = 0.502750, 4.5e-7 from pyNei's 0.5027494. The half called
genotype rule has no outside check either: plink2 refuses a VCF with a half
call unless `--vcf-half-call` tells it what to do with one, and the panel has
none.

Against pyNei: both libraries run `calc_per_var_distribs(stats="exp_het")` on
`tests/reference/sim_missing.vars` with the three populations,
`min_num_samples=20` and both values of `unbiased_exp_het`. The histogram counts
have to be equal and the mean of each population equal within 1e-12 relative,
because numpy and the Rust loop add the alleles of a variant and the variants
of a population in different orders.

The worked example, which becomes the first cargo test: 3 variants, 5 diploid
samples, pop1 = s1, s2 and pop2 = s3, s4, s5, `min_num_samples` 1.

| variant | genotypes | plain pop1 | plain pop2 | unbiased pop1 | unbiased pop2 |
|---|---|---|---|---|---|
| 1 | 0/0 2/1 0/0 0/0 0/. | 0.625 | 0 | 0.833333 | 0 |
| 2 | 0/0 0/0 0/1 1/0 ./. | 0 | 0.5 | 0 | 0.666667 |
| 3 | ./. ./. ./. ./. ./. | none | none | none | none |

Variant 1 in pop2 has c = 5, four alleles from two genotypes and one from the
half called one, all of allele 0, so both values are 0. The means over the
variants with a value are 0.3125 and 0.25 plain, 0.416667 and 0.333333
unbiased. The per variant values are the ones `test_calc_exp_het` asserts; the
means come from a run of pyNei at commit ef0ca6e with these two populations.

## The Rust interface

How the expected heterozygosity of a variant is worked out, and whether the
variant has one at all. The ploidies are `usize`, as they are in the
architecture, and they are two fields because that is what Open 3 recommends;
if the `ploidy` argument goes, they become one.

```rust
pub struct ExpHet {
    /// k, the exponent of the frequencies: the ploidy of the variants unless
    /// the caller asks for another one.
    pub exponent: usize,
    /// The ploidy of the variants, which turns the alleles a pop has called
    /// into genotypes for the min_num_samples test.
    pub ploidy: usize,
    /// How many called genotypes a pop needs at a variant to get a value.
    pub min_num_samples: u32,
    pub unbiased: bool,
}

impl ExpHet {
    /// The expected heterozygosity of one variant in one pop. `counts[a]` is
    /// how often allele a was called in the pop at this variant, and
    /// `called_alleles` is their sum; a missing allele is in neither. None in
    /// three cases: the pop has called fewer than `min_num_samples` genotypes;
    /// it has called nothing at all at this variant; or the unbiased one was
    /// asked for and the pop has called so little that 2n - 1 is 0, which for
    /// a diploid pop is one allele. A None counts for nothing: the variant is
    /// out of the mean and in no histogram bin.
    pub fn of_var(&self, counts: &[u32], called_alleles: u32) -> Option<f64>;
}
```

`counts` is what the row helper of the `variant` module fills for the sample
indices of one population, and the caller hands it the same buffer for every
variant. `of_var` compares `called_alleles` with `min_num_samples` times the
ploidy rather than dividing as pyNei does, which is the same test because the
ploidy is positive, and it never has to hold the 4.5 genotypes of a population
with a half called one.

## Open points

The owner decides these three. Until then the implementer follows the
"meanwhile" of each.

**Open 1: a chunk in which no genotype is called.** pyNei gives every variant
of it a plain expected heterozygosity of 1 and an unbiased one of -0.0, and the
histogram counts them, where the same variant among called ones gets no value.
The options are to reproduce it, which costs popnei a notion of the alleles of
a whole block that this calculation does not otherwise need and makes the value
of a variant depend on where the block boundaries fell, or to give no value
whenever the population has no called allele, which is what a record level
implementation does by itself and changes what a user sees only for a block
with nothing called in any sample. Recommendation: give no value. Meanwhile the
implementer writes it that way; the comparison against pyNei is safe, since the
reference panel drops 3 in 100 genotypes of 200 samples and no block of it is
without called genotypes.

**Open 2: the unbiased factor when the ploidy is not 2.** pyNei multiplies by
2n/(2n-1) at every ploidy. At ploidy 2 that is c/(c-1) and the estimator is
Nei's. At ploidy 4 it is neither: 1.1685 for the variant of three samples
above, a value over 1 that the mean counts and the histogram drops. The options
are to reproduce pyNei, which keeps every current result and leaves values above
1; to multiply by c/(c-1), which is the right correction for a pair of copies
but not for the k copies that 1 - sum p^k is about; or to refuse the unbiased
one above ploidy 2, which makes the default argument raise for every tetraploid
user. Recommendation: reproduce pyNei and say in the doc comment that the factor
is the diploid one, since no tetraploid result of pyNei is verified against
anything and nothing is built on this yet. Meanwhile the implementer reproduces
pyNei.

**Open 3: the `ploidy` argument.** When it differs from the data's ploidy,
pyNei uses it as the exponent and as the count of alleles the samples are
expected to hold, and not in the `min_num_samples` test. The options are to
reproduce all three uses, which in popnei takes deliberate work because the
count of called alleles comes from the row; to keep the argument as the
exponent alone and take the called alleles from the data, which is the `ExpHet`
above; or to drop the argument and always use the ploidy of the variants, which
loses nothing if nobody passes it. Recommendation: keep it as the exponent
alone. It is the same as pyNei whenever the argument matches the data, and it
is what the record level code gives for free. Meanwhile the implementer writes
that and compares against pyNei only at the data's own ploidy.

## Not in this spec

- `calc_per_var_distribs` itself, the mean, the histogram, the threads and how
  fast the per variant pass has to be: the item of this spec for that function,
  which the observed heterozygosity and the major allele frequency share.
- The allele counts and the frequencies per population that this calculation
  consumes, and what `min_num_samples` does to them: their own item of this
  spec.
- The observed heterozygosity, the polymorphism ratio and the per sample
  statistics: other items of this spec.
- The filter that drops variants by their observed heterozygosity:
  `docs/specs/filters.md`.
- The expected heterozygosities inside Jost's D, the differentiation between
  two populations: `_calc_pairwise_dest` of `pynei/dists.py` builds its own
  1 - sum p^k per population and over the pooled frequencies, and corrects both
  with the harmonic mean of the called genotypes of the two populations, which
  is not the correction here. `docs/specs/dists.md`.
- Fis, Fst and the other diversity summaries: pyNei does not have them and
  popnei does not add them.
