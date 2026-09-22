//! The statistics of the variants and of the individuals, per population.
//!
//! A population is a named set of individuals that a calculation treats as
//! a group, and every statistic of this module is calculated for each
//! population over its individuals alone. [`Pops`] is what a pass works
//! with: the name of each population and the indices of its individuals
//! among the individuals the pass gives, which are those of the source
//! after the filter of individuals of `docs/specs/filters.md` when the
//! variants carry one.
//!
//! `docs/specs/stats.md` has the design, and the row `stats` of section 9
//! of `docs/architecture.md` where the module sits.

use crate::block::{Block, BlockReader};
use crate::error::{Error, Result};
use crate::filters::resolve_individuals;
use crate::io::vcf::MAX_PLOIDY;
use crate::variant::{
    AlleleCounts, GtCounts, Needs, count_alleles, count_alleles_of, count_gts, count_gts_of,
};

/// The name of the one population of a calculation that was given no
/// populations, inherited from pyNei's `DEF_POP_NAME`.
pub const DEFAULT_POP_NAME: &str = "pop";

/// How many called genotypes a population needs at a variant for the
/// variant to have a value there when the caller asks for no other number,
/// 20, inherited from pyNei's `MIN_NUM_SAMPLES_FOR_POP_STAT`. Nobody has
/// measured whether 20 is the right threshold.
pub const DEFAULT_MIN_NUM_INDIVIDUALS: u32 = 20;

/// The major allele frequency below which a variant is polymorphic in a
/// population when the caller asks for no other number, 0.95, inherited
/// from pyNei's `DEF_POLY_THRESHOLD`. Nobody has measured whether 0.95 is
/// the right threshold.
pub const DEFAULT_POLY_THRESHOLD: f64 = 0.95;

/// The two ends of the range a histogram of a statistic covers when the
/// caller asks for no other, 0 and 1, which is where the five statistics
/// of a variant lie. It is pyNei's `default_range` of `_prepare_bins`.
pub const DEFAULT_HIST_RANGE: (f64, f64) = (0.0, 1.0);

/// How many bins a histogram of a statistic has when the caller asks for no
/// other number, 40, inherited from pyNei's `_prepare_bins`.
pub const DEFAULT_NUM_BINS: usize = 40;

/// The most bins a histogram of a statistic is built with, 100000.
///
/// A histogram a person reads has tens of bins, and [`DEFAULT_NUM_BINS`],
/// pyNei's, is 40. Every bin is a count of 8 bytes for each population and
/// each statistic, held once by a pass and once more by every chunk of rows
/// a thread is reading, so 100000 bins of the four statistics that have
/// them, over one population, are 3.2 MB in a chunk. More bins than this
/// are refused, because their counts are a vector no machine gives: the
/// allocation of 2^60 of them panics, and a panic in the core is a
/// `PanicException` in Python, which derives from `BaseException` and ends
/// a notebook, and the trap that ends the module in wasm.
pub const MAX_NUM_BINS: usize = 100_000;

/// The name a Python and a TypeScript user writes for bins of equal width,
/// which [`HistBins::linear`] builds. pyNei spells it `lineal`, the Spanish
/// word, and popnei refuses that name as any other unknown one, which the
/// owner decided on 22 September 2026.
pub const LINEAR_BINS: &str = "linear";

/// The name a Python and a TypeScript user writes for bins of equal ratio,
/// which [`HistBins::logarithmic`] builds.
pub const LOGARITHMIC_BINS: &str = "logarithmic";

/// The bins a histogram of a statistic has when the caller names no kind,
/// those of equal width, which cover their range evenly where the values of
/// a statistic lie.
pub const DEFAULT_BIN_TYPE: &str = LINEAR_BINS;

/// One population: its name and the indices of its individuals.
#[derive(Debug)]
struct Pop {
    name: String,
    /// The index of each individual among those of the reader, in the order
    /// the user named them.
    individuals: Vec<usize>,
    /// Whether those indices are every individual of the reader, in its
    /// order.
    is_all: bool,
}

/// The populations a pass calculates its statistics for: each one's name
/// and the indices of its individuals among the individuals of the reader.
///
/// The populations are in the order they were given, which is the order the
/// keys of a user's `pops` iterate in, and every result of this module
/// holds its values in that order. An individual can be in two
/// populations, and is in each of them once.
#[derive(Debug)]
pub struct Pops {
    pops: Vec<Pop>,
}

impl Pops {
    /// One population, named [`DEFAULT_POP_NAME`], of every individual of a
    /// reader of `num_individuals`, in the order of the reader: what a
    /// calculation works over when the user named no population.
    #[must_use]
    pub fn all(num_individuals: usize) -> Pops {
        Pops {
            pops: vec![Pop {
                name: DEFAULT_POP_NAME.to_owned(),
                individuals: (0..num_individuals).collect(),
                is_all: true,
            }],
        }
    }

    /// The populations a user named, each name looked up among
    /// `individuals`, the individuals the pass gives.
    ///
    /// `pops` is the name of each population with the names of its
    /// individuals, in the order the user gave them; it comes from a dict
    /// in Python and from an object in TypeScript, which hold each
    /// population name once.
    ///
    /// # Errors
    ///
    /// A name that is not one of `individuals`, a name that is twice in one
    /// population, a population that names no individual, and no population
    /// at all. The first three name the population, and the first two the
    /// name the user wrote.
    pub fn from_names(pops: &[(String, Vec<String>)], individuals: &[String]) -> Result<Pops> {
        if pops.is_empty() {
            return Err(Error::NoPop);
        }
        let mut of_the_pass = Vec::with_capacity(pops.len());
        for (name, named) in pops {
            let of_the_pop =
                resolve_individuals(named, individuals).map_err(|error| in_the_pop(error, name))?;
            of_the_pass.push(Pop {
                name: name.clone(),
                is_all: every_individual_in_order(&of_the_pop, individuals.len()),
                individuals: of_the_pop,
            });
        }
        Ok(Pops { pops: of_the_pass })
    }

    /// How many populations there are, which is 1 when the user named none.
    #[must_use]
    pub fn len(&self) -> usize {
        self.pops.len()
    }

    /// Whether there is no population, which the constructors do not build:
    /// `pops` with no population is refused, and [`Pops::all`] gives one.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.pops.is_empty()
    }

    /// The name of one population, as the user wrote it.
    ///
    /// `pop` is a population of `0..len()`, which is how a caller walks
    /// them; a number at or beyond `len()` is no population of this and has
    /// no name here.
    #[must_use]
    pub fn name(&self, pop: usize) -> &str {
        self.pops.get(pop).map_or("", |pop| pop.name.as_str())
    }

    /// The index of each individual of one population among the individuals
    /// of the reader, in the order the user named them.
    ///
    /// `pop` is a population of `0..len()`; a number at or beyond `len()`
    /// is no population of this and has no individual here.
    #[must_use]
    pub fn individuals(&self, pop: usize) -> &[usize] {
        self.pops
            .get(pop)
            .map_or(&[][..], |pop| pop.individuals.as_slice())
    }

    /// Whether the population is every individual of the reader in the
    /// order of the reader, which a caller that reads a row of genotypes as
    /// it is asks before it counts.
    ///
    /// `pop` is a population of `0..len()`; a number at or beyond `len()`
    /// is no population of this and is not every individual.
    #[must_use]
    pub fn is_all(&self, pop: usize) -> bool {
        self.pops.get(pop).is_some_and(|pop| pop.is_all)
    }
}

/// The error of [`resolve_individuals`] over the names of one population,
/// with the population they were given in.
///
/// The filter of individuals and a population are given the same three
/// refusals by the same function, and what a user has to look at differs:
/// there the names are the whole of what they wrote, and here they are the
/// names of one of their populations, which is the one the message has to
/// name. An error of another kind, which that function does not give,
/// travels on as it is.
#[expect(
    clippy::wildcard_enum_match_arm,
    reason = "`resolve_individuals` fails with these three cases alone, and an error of \
              any other case of the crate is left as it is rather than named after a \
              population it may have nothing to do with"
)]
fn in_the_pop(error: Error, pop: &str) -> Error {
    match error {
        Error::IndividualNotInTheSource { name } => Error::IndividualOfAPopNotInThePass {
            pop: pop.to_owned(),
            name,
        },
        Error::IndividualNamedTwice { name } => Error::IndividualNamedTwiceInAPop {
            pop: pop.to_owned(),
            name,
        },
        Error::NoIndividualNamed => Error::PopWithNoIndividual {
            pop: pop.to_owned(),
        },
        of_another_kind => of_another_kind,
    }
}

/// Whether `individuals` is every individual of a reader of
/// `num_individuals`, in the order of the reader: the population a caller
/// counts by reading a row of genotypes as it is.
fn every_individual_in_order(individuals: &[usize], num_individuals: usize) -> bool {
    individuals.len() == num_individuals
        && individuals
            .iter()
            .enumerate()
            .all(|(of_the_reader, individual)| of_the_reader == *individual)
}

/// The edges of the bins that a histogram counts the values of a statistic
/// in.
///
/// There is one edge more than there are bins. A value falls in the bin
/// whose left edge is at most the value and whose right edge is above it,
/// and the last bin takes its right edge too, as `numpy.histogram` does, so
/// an observed heterozygosity of exactly 1 is in the last bin. A value
/// outside the range is in no bin.
///
/// The edges are the ones numpy's `linspace` computes, the start of the
/// range plus i times the width of a bin for the i-th of them and the end
/// of the range for the last, so that a value that falls on an edge falls
/// on the same side of it in popnei and in pyNei. The edges of bins of
/// equal ratio are 10 raised to the same over the base 10 logarithms of the
/// two ends, which is what numpy's `logspace` computes.
#[derive(Debug, Clone)]
pub struct HistBins {
    /// The edges, from the start of the range up to its end, one more than
    /// there are bins.
    edges: Vec<f64>,
}

impl HistBins {
    /// `num_bins` bins of equal width from `start` to `end`.
    ///
    /// # Errors
    ///
    /// When `num_bins` is 0 or above [`MAX_NUM_BINS`], when `start` is not
    /// below `end`, when either of the two is NaN or infinite, which leaves
    /// every edge between them NaN, and when the two are further apart than
    /// a float64 goes, which leaves the width of a bin infinite.
    pub fn linear(start: f64, end: f64, num_bins: usize) -> Result<HistBins> {
        check_the_range(start, end, num_bins)?;
        let width_of_a_bin = (end - start) / num_bins as f64;
        Ok(HistBins {
            edges: edges_of(start, width_of_a_bin, num_bins, end),
        })
    }

    /// `num_bins` bins of equal ratio from `start` to `end`: each edge is
    /// the one before it times a fixed factor.
    ///
    /// # Errors
    ///
    /// Those of [`HistBins::linear`], and a `start` of 0 or below, which no
    /// factor takes anywhere.
    pub fn logarithmic(start: f64, end: f64, num_bins: usize) -> Result<HistBins> {
        check_the_range(start, end, num_bins)?;
        if start <= 0.0 {
            return Err(Error::HistLogRangeNotAboveZero { start });
        }
        // numpy's `logspace` raises 10 to the edges of equal width between
        // the two logarithms, the last of which is the logarithm of the end
        // and not the end itself.
        let (first, last) = (start.log10(), end.log10());
        let width_of_a_bin = (last - first) / num_bins as f64;
        let mut edges = edges_of(first, width_of_a_bin, num_bins, last);
        for edge in &mut edges {
            *edge = 10_f64.powf(*edge);
        }
        Ok(HistBins { edges })
    }

    /// The edges, from the start of the range up to its end, one more than
    /// there are bins.
    #[must_use]
    pub fn edges(&self) -> &[f64] {
        &self.edges
    }

    /// How many bins there are, which is one less than the edges.
    #[must_use]
    pub fn num_bins(&self) -> usize {
        // A `HistBins` holds 2 edges at least: both constructors refuse a
        // histogram of no bin.
        self.edges.len().saturating_sub(1)
    }

    /// The bin `value` falls in, or `None` when it is outside the range and
    /// in no bin.
    ///
    /// A value falls in the bin whose left edge is at most the value and
    /// whose right edge is above it, and the last bin takes its right edge
    /// too, which is what `numpy.histogram` does.
    #[must_use]
    pub fn bin_of(&self, value: f64) -> Option<usize> {
        let first = *self.edges.first()?;
        let last = *self.edges.last()?;
        if value.is_nan() || value < first || value > last {
            return None;
        }
        // The value is at the first edge or above it, so one edge at least
        // is at most the value and the count below is 1 or more: the bin
        // that starts at the last of those edges is the one the value falls
        // in. A value at the end of the range has every edge at or below it
        // and falls in the last bin, which takes its right edge.
        let edges_at_or_below = self.edges.partition_point(|edge| *edge <= value);
        Some(
            edges_at_or_below
                .saturating_sub(1)
                .min(self.num_bins().saturating_sub(1)),
        )
    }
}

/// The edges of `num_bins` bins over a range: the start of the range plus
/// i times the width of a bin for the i-th of them, and the end of the
/// range for the last, which is what numpy's `linspace` computes.
///
/// The last edge is the end itself and not the start plus `num_bins` times
/// the width, so that the range the bins cover is the one that was asked
/// for whatever the rounding of the widths added up to.
fn edges_of(start: f64, width_of_a_bin: f64, num_bins: usize, end: f64) -> Vec<f64> {
    // One edge more than there are bins; a `num_bins` that saturates asks
    // for more memory than a machine gives either way.
    let mut edges = Vec::with_capacity(num_bins.saturating_add(1));
    edges.extend((0..num_bins).map(|bin| start + bin as f64 * width_of_a_bin));
    edges.push(end);
    edges
}

/// The checks both constructors of [`HistBins`] make on the range and the
/// number of bins they were given.
///
/// # Errors
///
/// A `num_bins` of 0 or above [`MAX_NUM_BINS`], a `start` that is not below
/// `end`, a `start` or an `end` that is NaN or infinite, and two ends
/// further apart than a float64 goes.
fn check_the_range(start: f64, end: f64, num_bins: usize) -> Result<()> {
    if num_bins == 0 {
        return Err(Error::HistWithNoBin);
    }
    if num_bins > MAX_NUM_BINS {
        return Err(Error::HistTooManyBins {
            num_bins,
            largest: MAX_NUM_BINS,
        });
    }
    if !start.is_finite() || !end.is_finite() || start >= end {
        return Err(Error::HistRangeNotGoingUp { start, end });
    }
    // The width of a bin is the distance between the two ends over the
    // bins, and two finite ends can be further apart than a float64 goes:
    // -1e308 to 1e308 leaves the width infinite and the edges NaN,
    // infinite, infinite, infinite and 1e308, which do not go up, so the
    // search for the bin of a value puts every value in the first bin.
    if !(end - start).is_finite() {
        return Err(Error::HistRangeTooWide { start, end });
    }
    Ok(())
}

/// The observed heterozygosity of one variant in one population: the share
/// of its called genotypes that are heterozygous.
///
/// A genotype is heterozygous when it is called and its alleles are not all
/// the same, at any ploidy, and a half called genotype is missing and
/// counts in neither the heterozygous ones nor the called ones. At ploidy 1
/// no genotype is heterozygous and every variant with a called genotype has
/// 0.
///
/// It mirrors `_calc_obs_het_per_var` of pyNei, which holds it to no
/// threshold of how much data a population has; popnei holds it to
/// `min_num_individuals` like the other statistics, which the owner decided
/// on 22 September 2026.
#[derive(Debug, Clone, Copy)]
pub struct ObsHet {
    /// How many genotypes the population has to have called at the variant
    /// for the variant to have a value there.
    min_num_individuals: u32,
}

impl ObsHet {
    /// The statistic, with how many called genotypes a population needs at
    /// a variant to have a value there.
    #[must_use]
    pub fn new(min_num_individuals: u32) -> ObsHet {
        ObsHet {
            min_num_individuals,
        }
    }

    /// The heterozygous genotypes over the called ones, from the counts of
    /// the genotypes of the population at one variant.
    ///
    /// `None` when the population has called nothing at the variant, and
    /// when it has called fewer genotypes than `min_num_individuals`: the
    /// variant is then out of the mean and in no bin of the histogram.
    #[must_use]
    pub fn of_var(&self, counts: GtCounts) -> Option<f64> {
        if counts.called == 0 || counts.called < self.min_num_individuals {
            return None;
        }
        Some(f64::from(counts.het) / f64::from(counts.called))
    }
}

/// The major allele frequency of one variant in one population: how often
/// its commonest allele was called among the called alleles of the
/// population.
///
/// "maf" is, in pyNei and in popnei, the frequency of the major allele,
/// where most of the literature and plink2 use the same letters for the
/// minor one. An allele is counted each time it was called, in a half
/// called genotype too, and every allele of a multiallelic variant has its
/// own count; when two alleles tie for the largest count the value is the
/// same whichever of them is called the major one. It mirrors
/// `_calc_maf_per_var` of pyNei.
#[derive(Debug, Clone, Copy)]
pub struct Maf {
    /// How many alleles the population has to have called at the variant,
    /// which is `min_num_individuals` genotypes of the ploidy.
    min_called_alleles: u64,
}

impl Maf {
    /// The statistic, with the ploidy of the variants and how many called
    /// genotypes a population needs at a variant to have a value there.
    ///
    /// # Errors
    ///
    /// A `ploidy` of 0 or above the largest ploidy a reader of popnei
    /// gives, 255.
    pub fn new(ploidy: usize, min_num_individuals: u32) -> Result<Maf> {
        let ploidy = checked_ploidy("ploidy", ploidy)?;
        Ok(Maf {
            min_called_alleles: min_called_alleles(min_num_individuals, ploidy),
        })
    }

    /// The largest of `counts` over `called_alleles`, which is their sum.
    ///
    /// `counts[a]` is how often the allele a was called in the population at
    /// this variant, which [`count_alleles_of`] fills. `None` when the
    /// population has called nothing at the variant, and when it has called
    /// fewer than `min_num_individuals` genotypes, which is `called_alleles`
    /// below that many times the ploidy.
    ///
    /// [`count_alleles_of`]: crate::variant::count_alleles_of
    #[must_use]
    pub fn of_var(&self, counts: &AlleleCounts, called_alleles: u32) -> Option<f64> {
        if called_alleles == 0 || u64::from(called_alleles) < self.min_called_alleles {
            return None;
        }
        let of_the_major_allele = counts.iter().copied().max()?;
        Some(f64::from(of_the_major_allele) / f64::from(called_alleles))
    }
}

/// The expected heterozygosity of one variant in one population, plain or
/// unbiased: the chance that gene copies taken at random from the
/// population are not all of the same allele.
///
/// The plain one is `1 - sum over a of p_a^k`, where p_a is the frequency of
/// the allele a among the alleles the population called at the variant and k
/// is the exponent. The unbiased one corrects for the frequencies being
/// estimated from the same copies the statistic is computed over, which
/// makes the plain one too small on average: with c the called alleles and
/// c_a the count of the allele a among them, it is
/// `1 - sum over a of (c_a (c_a - 1) ... (c_a - k + 1)) / (c (c - 1) ... (c - k + 1))`,
/// the chance that k copies drawn from the called ones without replacement
/// are all alike, taken from 1. At k = 2 that is Nei's 1978 correction,
/// `(c / (c - 1))` times the plain one, which GenAlEx prints as the
/// unbiased heterozygosity for codominant data.
///
/// It mirrors `_calc_exp_het_per_var` and `_calc_unbiased_exp_het_per_var`
/// of pyNei, which apply the factor of a diploid population at every
/// ploidy; popnei applies the correction of the exponent in hand, which the
/// owner decided on 22 September 2026.
#[derive(Debug, Clone, Copy)]
pub struct ExpHet {
    /// k, the exponent the allele frequencies are raised to and the number
    /// of factors of the products of the unbiased one, 1 to 255.
    exponent: u32,
    /// How many alleles the population has to have called at the variant,
    /// which is `min_num_individuals` genotypes of the ploidy of the
    /// variants. The exponent is never used for it.
    min_called_alleles: u64,
}

impl ExpHet {
    /// The statistic, with the exponent, the ploidy of the variants and how
    /// many called genotypes a population needs at a variant to have a
    /// value there.
    ///
    /// `exponent` is k, the number the allele frequencies are raised to and
    /// the number of factors of the products of the unbiased one: the ploidy
    /// of the variants unless the caller asks for another one. `ploidy` is
    /// the ploidy of the variants, which turns the alleles a population has
    /// called into called genotypes for the `min_num_individuals` test; the
    /// exponent is never used for that.
    ///
    /// # Errors
    ///
    /// An `exponent` or a `ploidy` of 0 or above the largest ploidy a
    /// reader of popnei gives, 255.
    pub fn new(exponent: usize, ploidy: usize, min_num_individuals: u32) -> Result<ExpHet> {
        let exponent = checked_ploidy("exponent", exponent)?;
        let ploidy = checked_ploidy("ploidy", ploidy)?;
        Ok(ExpHet {
            exponent,
            min_called_alleles: min_called_alleles(min_num_individuals, ploidy),
        })
    }

    /// The expected heterozygosity of one variant in one population, the
    /// plain one or, with `unbiased`, the unbiased one.
    ///
    /// `counts[a]` is how often the allele a was called in the population at
    /// this variant, which [`count_alleles_of`] fills, and `called_alleles`
    /// is their sum; a missing allele is in neither. `None` in three cases:
    /// the population has called fewer than `min_num_individuals` genotypes;
    /// it has called nothing at all at this variant; or the unbiased one was
    /// asked for and `called_alleles` is below the exponent, which for a
    /// diploid population is one called allele.
    ///
    /// [`count_alleles_of`]: crate::variant::count_alleles_of
    #[must_use]
    pub fn of_var(
        &self,
        counts: &AlleleCounts,
        called_alleles: u32,
        unbiased: bool,
    ) -> Option<f64> {
        if called_alleles == 0 || u64::from(called_alleles) < self.min_called_alleles {
            return None;
        }
        if unbiased {
            if called_alleles < self.exponent {
                return None;
            }
            return Some(1.0 - self.all_alike_without_replacement(counts, called_alleles));
        }
        Some(1.0 - self.all_alike_with_replacement(counts, called_alleles))
    }

    /// The chance that k gene copies drawn from the called ones with
    /// replacement are all of the same allele, the sum over the alleles of
    /// `p_a^k`, which the plain one is taken from 1.
    ///
    /// The alleles are added in the order of their number, so that two runs
    /// over the same variant give the same bits.
    fn all_alike_with_replacement(&self, counts: &AlleleCounts, called_alleles: u32) -> f64 {
        let called_alleles = f64::from(called_alleles);
        counts
            .iter()
            .filter(|count| **count > 0)
            .map(|count| raised(f64::from(*count) / called_alleles, self.exponent))
            .sum()
    }

    /// The chance that k gene copies drawn from the called ones without
    /// replacement are all of the same allele, the sum over the alleles of
    /// `(c_a (c_a - 1) ... (c_a - k + 1)) / (c (c - 1) ... (c - k + 1))`,
    /// which the unbiased one is taken from 1.
    ///
    /// Each term is the product of its k factors one over another and not
    /// one product over the other, because `c (c - 1) ... (c - k + 1)` of a
    /// population of a million called alleles and an exponent of 255 is
    /// above the largest float64 and would leave the term NaN. Every factor
    /// is a number from 0 to 1, and an allele called fewer than k times has
    /// a factor of 0 and a term of 0: there are not k copies of it to draw.
    /// The caller has checked that `called_alleles` is k or more, so no
    /// denominator is 0.
    fn all_alike_without_replacement(&self, counts: &AlleleCounts, called_alleles: u32) -> f64 {
        let called_alleles = f64::from(called_alleles);
        counts
            .iter()
            .filter(|count| **count >= self.exponent)
            .map(|count| {
                let count = f64::from(*count);
                let mut all_alike = 1.0;
                let mut drawn = 0.0;
                for _ in 0..self.exponent {
                    all_alike *= (count - drawn) / (called_alleles - drawn);
                    drawn += 1.0;
                }
                all_alike
            })
            .sum()
    }
}

/// `value` raised to the power `times`, as `times` multiplications.
///
/// `powi` and `powf` are not rounded the same on every platform, and a
/// multiplication is, so the plain expected heterozygosity of a diploid
/// population is the `p * p` numpy computes, to the bit, on every machine
/// popnei runs on. `times` is 255 at most, the largest exponent a statistic
/// is built with.
fn raised(value: f64, times: u32) -> f64 {
    let mut raised = 1.0;
    for _ in 0..times {
        raised *= value;
    }
    raised
}

/// `value` as the ploidy or the exponent of a statistic of one variant: 1
/// at least and [`MAX_PLOIDY`] at most, the largest ploidy a reader of
/// popnei gives.
///
/// `kind` is which of the two it is, `ploidy` or `exponent`, which the
/// error names.
///
/// # Errors
///
/// A `value` of 0 or above [`MAX_PLOIDY`].
fn checked_ploidy(kind: &'static str, value: usize) -> Result<u32> {
    let out_of_range = || Error::StatPloidyOutOfRange {
        kind,
        value,
        largest: MAX_PLOIDY,
    };
    if value == 0 || value > MAX_PLOIDY {
        return Err(out_of_range());
    }
    u32::try_from(value).map_err(|_| out_of_range())
}

/// How many alleles a population has to have called at a variant to have a
/// value there: `min_num_individuals` genotypes of `ploidy` alleles each.
///
/// The three statistics compare the called alleles with this number rather
/// than dividing them by the ploidy as pyNei does, which is the same test
/// because the ploidy is positive, and which never has to hold the 4.5
/// genotypes of a population with a half called one.
fn min_called_alleles(min_num_individuals: u32, ploidy: u32) -> u64 {
    // At most 4295 million genotypes of a ploidy of 255 at most, which is
    // 1.1e12 and fits in a u64 many times over; the saturation is there
    // because the plain operator does not compile, and a threshold that
    // saturated would be one no population ever meets, which is what a
    // number that large asks for.
    u64::from(min_num_individuals).saturating_mul(u64::from(ploidy))
}

/// One of the five statistics that [`calc_per_var_distribs`] calculates for
/// every variant and every population.
///
/// The names are the ones a Python and a TypeScript user writes, and each
/// is the field of the result that holds its distribution or its counts.
/// The expected heterozygosity, plain, and the unbiased one are two
/// statistics and not one with a switch, which the owner decided on 22
/// September 2026, so that a user asks for the two like any two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerVarStat {
    /// The heterozygous genotypes of a population over its called ones.
    ObsHet,
    /// The commonest allele of a population over its called alleles.
    Maf,
    /// The chance that gene copies taken at random from a population are
    /// not all of the same allele, from the frequencies as they are.
    ExpHet,
    /// The same corrected for the frequencies being estimated from the
    /// copies the statistic is computed over.
    UnbiasedExpHet,
    /// How many of the variants vary in a population, in three counts and
    /// two ratios, which is a count and not a distribution.
    PolyVarsRatio,
}

/// What one pass of [`calc_per_var_distribs`] calculates: which statistics,
/// for which populations, in which bins and with the thresholds of each.
#[derive(Debug)]
pub struct PerVarDistribsConfig {
    /// The statistics to calculate. A statistic that is named twice is
    /// calculated once, and one that is named by nobody is `None` in the
    /// result; asking for fewer is a saving of work and changes no value.
    pub stats: Vec<PerVarStat>,
    /// The populations, each of which gets its own value of every
    /// statistic, in the order they were given.
    pub pops: Pops,
    /// The bins the values of every distribution are counted in.
    pub bins: HistBins,
    /// The observed heterozygosity of one variant in one population.
    pub obs_het: ObsHet,
    /// The major allele frequency, which the polymorphism ratio counts
    /// too.
    pub maf: Maf,
    /// The expected heterozygosity, which gives the plain value and the
    /// unbiased one.
    pub exp_het: ExpHet,
    /// Below this major allele frequency a variant is polymorphic in a
    /// population. A number from 0 to 1, both included; anything else is
    /// an error of [`calc_per_var_distribs`].
    pub poly_threshold: f64,
}

/// The distribution of one statistic over the variants of a pass, for each
/// population: the mean of the variants that had a value and how many of
/// them fell in each bin.
///
/// A variant that has no value in a population, because the population has
/// too little data at it, is out of the mean and in no bin, so the
/// histograms of two populations can count different numbers of variants. A
/// value outside the range of the bins is in the mean and in no bin, which
/// is what `numpy.histogram` leaves out too.
#[derive(Debug, Clone)]
pub struct StatsDistrib {
    /// The bins the values were counted in, the same for every population.
    bins: HistBins,
    /// What the pass added up for each population, in the order of
    /// [`Pops`].
    pops: Vec<Accumulated>,
}

impl StatsDistrib {
    /// How many populations the distribution holds a value for, which is
    /// the populations of the pass.
    #[must_use]
    pub fn num_pops(&self) -> usize {
        self.pops.len()
    }

    /// The mean of the variants of one population that had a value, and
    /// `None` when none of them had one.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no mean here.
    #[must_use]
    pub fn mean(&self, pop: usize) -> Option<f64> {
        let of_the_pop = self.pops.get(pop)?;
        if of_the_pop.num_vars_with_value == 0 {
            return None;
        }
        // Every count of popnei is below 2^53, where a `f64` holds the
        // whole numbers exactly, so the division is the sum over the count.
        Some(of_the_pop.sum / of_the_pop.num_vars_with_value as f64)
    }

    /// How many variants of one population had a value: those the mean is
    /// over.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and had no variant here.
    #[must_use]
    pub fn num_vars_with_value(&self, pop: usize) -> u64 {
        self.pops
            .get(pop)
            .map_or(0, |of_the_pop| of_the_pop.num_vars_with_value)
    }

    /// How many variants of one population fell in each bin, one count for
    /// each bin of [`StatsDistrib::bins`].
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no counts here.
    #[must_use]
    pub fn hist_counts(&self, pop: usize) -> &[u64] {
        self.pops
            .get(pop)
            .map_or(&[][..], |of_the_pop| of_the_pop.hist.as_slice())
    }

    /// The bins the values were counted in, whose edges every result of a
    /// pass gives back.
    #[must_use]
    pub fn bins(&self) -> &HistBins {
        &self.bins
    }
}

/// How many of the variants of a pass vary in each population: the three
/// counts of the polymorphism ratio and the two ratios that follow from
/// them.
///
/// A variant is polymorphic in a population when its major allele frequency
/// there is below the threshold of the pass, strictly, and variable when
/// that frequency is below 1. Both are counted among the variants that have
/// a major allele frequency in the population, those at which it called
/// `min_num_individuals` genotypes or more.
#[derive(Debug, Clone)]
pub struct PolyVarsStats {
    /// The three counts of each population, in the order of [`Pops`].
    pops: Vec<PolyCounts>,
}

impl PolyVarsStats {
    /// How many populations it holds the counts of.
    #[must_use]
    pub fn num_pops(&self) -> usize {
        self.pops.len()
    }

    /// The variants that are polymorphic in one population.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no count here.
    #[must_use]
    pub fn num_poly(&self, pop: usize) -> u64 {
        self.pops.get(pop).map_or(0, |counts| counts.num_poly)
    }

    /// The variants that are variable in one population.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no count here.
    #[must_use]
    pub fn num_variable(&self, pop: usize) -> u64 {
        self.pops.get(pop).map_or(0, |counts| counts.num_variable)
    }

    /// The variants that have a major allele frequency in one population,
    /// which the other two counts are among.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no count here.
    #[must_use]
    pub fn num_vars_with_data(&self, pop: usize) -> u64 {
        self.pops
            .get(pop)
            .map_or(0, |counts| counts.num_vars_with_data)
    }

    /// The polymorphic variants over the ones with data, and `None` when
    /// the population has none with data.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no ratio here.
    #[must_use]
    pub fn poly_ratio(&self, pop: usize) -> Option<f64> {
        let counts = self.pops.get(pop)?;
        ratio_of(counts.num_poly, counts.num_vars_with_data)
    }

    /// The polymorphic variants over the variable ones, and `None` when the
    /// population has no variable one.
    ///
    /// `pop` is a population of `0..num_pops()`; a number at or beyond
    /// `num_pops()` is no population of this and has no ratio here.
    #[must_use]
    pub fn poly_ratio_over_variables(&self, pop: usize) -> Option<f64> {
        let counts = self.pops.get(pop)?;
        ratio_of(counts.num_poly, counts.num_variable)
    }
}

/// `of_them` over `of_all`, and `None` when `of_all` is 0, which every
/// ratio of the polymorphism counts is then.
fn ratio_of(of_them: u64, of_all: u64) -> Option<f64> {
    if of_all == 0 {
        return None;
    }
    // Every count of popnei is below 2^53, where a `f64` holds the whole
    // numbers exactly.
    Some(of_them as f64 / of_all as f64)
}

/// What one pass of [`calc_per_var_distribs`] gives back: the distribution
/// of each statistic that was asked for, per population, and how many
/// variants the pass gave.
///
/// A statistic that was not asked for is `None`. The populations of every
/// field are in the order of the [`Pops`] of the pass, which is the order
/// the user named them in.
#[derive(Debug)]
pub struct PerVarDistribs {
    /// The observed heterozygosity.
    pub obs_het: Option<StatsDistrib>,
    /// The major allele frequency.
    pub maf: Option<StatsDistrib>,
    /// The plain expected heterozygosity.
    pub exp_het: Option<StatsDistrib>,
    /// The unbiased expected heterozygosity.
    pub unbiased_exp_het: Option<StatsDistrib>,
    /// The counts of the polymorphism ratio.
    pub poly_vars_ratio: Option<PolyVarsStats>,
    /// The variants the pass gave, after the steps the variants carried.
    pub num_vars: u64,
}

/// How many rows of a block one chunk of the pass reads.
///
/// The rows of a block are added up chunk by chunk and the chunks are added
/// together in the order of the block, so the sum of a statistic does not
/// depend on how many threads read the block, which rayon's own `sum` would
/// make it: it joins the parts in an order it chooses at run time. The
/// number is fixed for the same reason, and 64 rows of 1000 diploid
/// individuals are 128000 genotypes, enough work for one task of rayon.
const ROWS_PER_CHUNK: usize = 64;

/// Which of the five statistics a pass was asked for, and what has to be
/// counted for them.
#[derive(Debug, Clone, Copy)]
struct Asked {
    obs_het: bool,
    maf: bool,
    exp_het: bool,
    unbiased_exp_het: bool,
    poly_vars_ratio: bool,
}

impl Asked {
    /// The statistics of a configuration, with a statistic named twice
    /// taken once.
    fn of(stats: &[PerVarStat]) -> Asked {
        let mut asked = Asked {
            obs_het: false,
            maf: false,
            exp_het: false,
            unbiased_exp_het: false,
            poly_vars_ratio: false,
        };
        for stat in stats {
            match *stat {
                PerVarStat::ObsHet => asked.obs_het = true,
                PerVarStat::Maf => asked.maf = true,
                PerVarStat::ExpHet => asked.exp_het = true,
                PerVarStat::UnbiasedExpHet => asked.unbiased_exp_het = true,
                PerVarStat::PolyVarsRatio => asked.poly_vars_ratio = true,
            }
        }
        asked
    }

    /// Whether the alleles of a row have to be counted for a population:
    /// three of the five statistics follow from those counts, and the
    /// fourth, the polymorphism ratio, from the major allele frequency they
    /// give.
    fn the_allele_counts(self) -> bool {
        self.maf || self.exp_het || self.unbiased_exp_het || self.poly_vars_ratio
    }

    /// Whether the major allele frequency of a row has to be worked out for
    /// a population: for itself, and for the polymorphism ratio, which
    /// counts the variants whose frequency is below a threshold.
    fn the_maf(self) -> bool {
        self.maf || self.poly_vars_ratio
    }
}

/// What a pass has added up of one statistic over one population: the sum
/// of the values, how many variants had one, and how many of them fell in
/// each bin.
#[derive(Debug, Clone)]
struct Accumulated {
    sum: f64,
    num_vars_with_value: u64,
    /// One count for each bin, and no count at all for a statistic that
    /// was not asked for, which nothing is added to.
    hist: Vec<u64>,
}

impl Accumulated {
    /// The accumulator of a statistic before any variant is read, with the
    /// bins of a statistic that was asked for and none of one that was not.
    fn of(num_bins: usize, asked_for: bool) -> Accumulated {
        Accumulated {
            sum: 0.0,
            num_vars_with_value: 0,
            hist: vec![0; if asked_for { num_bins } else { 0 }],
        }
    }

    /// It adds the value one variant had in this population.
    fn add(&mut self, value: f64, bins: &HistBins) {
        self.sum += value;
        // One variant of the pass, and a pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds.
        self.num_vars_with_value = self.num_vars_with_value.saturating_add(1);
        if let Some(bin) = bins.bin_of(value)
            && let Some(count) = self.hist.get_mut(bin)
        {
            *count = count.saturating_add(1);
        }
    }

    /// It empties the accumulator, so that one chunk of rows after another
    /// is added up in it without its bins being allocated again.
    fn forget_what_it_holds(&mut self) {
        self.sum = 0.0;
        self.num_vars_with_value = 0;
        self.hist.fill(0);
    }

    /// It adds what one chunk of rows found to what the pass has.
    ///
    /// The chunks are added in the order of the block, so the sum does not
    /// depend on how many threads read them.
    fn add_the_chunk(&mut self, of_the_chunk: &Accumulated) {
        self.sum += of_the_chunk.sum;
        self.num_vars_with_value = self
            .num_vars_with_value
            .saturating_add(of_the_chunk.num_vars_with_value);
        for (count, of_the_chunk) in self.hist.iter_mut().zip(&of_the_chunk.hist) {
            *count = count.saturating_add(*of_the_chunk);
        }
    }
}

/// The three counts of the polymorphism ratio of one population.
#[derive(Debug, Clone, Copy)]
struct PolyCounts {
    num_poly: u64,
    num_variable: u64,
    num_vars_with_data: u64,
}

impl PolyCounts {
    /// The counts before any variant is read.
    fn none() -> PolyCounts {
        PolyCounts {
            num_poly: 0,
            num_variable: 0,
            num_vars_with_data: 0,
        }
    }

    /// It counts one variant that has a major allele frequency in this
    /// population, polymorphic when that frequency is below `threshold` and
    /// variable when it is below 1.
    fn add(&mut self, maf: f64, threshold: f64) {
        // One variant of the pass each, which no source holds
        // 18446744073709551615 of.
        self.num_vars_with_data = self.num_vars_with_data.saturating_add(1);
        if maf < 1.0 {
            self.num_variable = self.num_variable.saturating_add(1);
        }
        if maf < threshold {
            self.num_poly = self.num_poly.saturating_add(1);
        }
    }

    /// It adds what one chunk of rows found to what the pass has.
    fn add_the_chunk(&mut self, of_the_chunk: PolyCounts) {
        self.num_poly = self.num_poly.saturating_add(of_the_chunk.num_poly);
        self.num_variable = self.num_variable.saturating_add(of_the_chunk.num_variable);
        self.num_vars_with_data = self
            .num_vars_with_data
            .saturating_add(of_the_chunk.num_vars_with_data);
    }
}

/// What a pass, or one chunk of the rows of a block, has added up of every
/// statistic over one population.
#[derive(Debug, Clone)]
struct OfAPop {
    obs_het: Accumulated,
    maf: Accumulated,
    exp_het: Accumulated,
    unbiased_exp_het: Accumulated,
    poly: PolyCounts,
}

/// What a pass, or one chunk of the rows of a block, has added up over
/// every population.
///
/// Its size is the populations times the bins of each statistic, and grows
/// neither with the variants of the pass nor with the individuals: it is
/// what is kept from one block to the next.
#[derive(Debug, Clone)]
struct Totals {
    pops: Vec<OfAPop>,
}

impl Totals {
    /// The accumulators of `num_pops` populations before any variant is
    /// read, with the bins of the statistics that were asked for.
    fn of(num_pops: usize, num_bins: usize, asked: Asked) -> Totals {
        Totals {
            pops: (0..num_pops)
                .map(|_| OfAPop {
                    obs_het: Accumulated::of(num_bins, asked.obs_het),
                    maf: Accumulated::of(num_bins, asked.maf),
                    exp_het: Accumulated::of(num_bins, asked.exp_het),
                    unbiased_exp_het: Accumulated::of(num_bins, asked.unbiased_exp_het),
                    poly: PolyCounts::none(),
                })
                .collect(),
        }
    }

    /// It adds what one chunk of rows found to what the pass has, one
    /// population at a time.
    fn add_the_chunk(&mut self, of_the_chunk: &Totals) {
        for (of_the_pass, of_the_chunk) in self.pops.iter_mut().zip(&of_the_chunk.pops) {
            of_the_pass.obs_het.add_the_chunk(&of_the_chunk.obs_het);
            of_the_pass.maf.add_the_chunk(&of_the_chunk.maf);
            of_the_pass.exp_het.add_the_chunk(&of_the_chunk.exp_het);
            of_the_pass
                .unbiased_exp_het
                .add_the_chunk(&of_the_chunk.unbiased_exp_het);
            of_the_pass.poly.add_the_chunk(of_the_chunk.poly);
        }
    }

    /// It empties every accumulator, so that one chunk of rows after
    /// another is added up in the same `Totals`.
    fn forget_what_it_holds(&mut self) {
        for of_the_pop in &mut self.pops {
            of_the_pop.obs_het.forget_what_it_holds();
            of_the_pop.maf.forget_what_it_holds();
            of_the_pop.exp_het.forget_what_it_holds();
            of_the_pop.unbiased_exp_het.forget_what_it_holds();
            of_the_pop.poly = PolyCounts::none();
        }
    }

    /// The distribution of one statistic, with the accumulator of each
    /// population that `of` picks out of it.
    fn distrib_of(&self, bins: &HistBins, of: impl Fn(&OfAPop) -> &Accumulated) -> StatsDistrib {
        StatsDistrib {
            bins: bins.clone(),
            pops: self.pops.iter().map(|pop| of(pop).clone()).collect(),
        }
    }

    /// The counts of the polymorphism ratio of every population.
    fn poly_vars_stats(&self) -> PolyVarsStats {
        PolyVarsStats {
            pops: self.pops.iter().map(|pop| pop.poly).collect(),
        }
    }
}

/// The distributions of the statistics of `config` over the variants
/// `reader` gives, which is one pass over the source through the steps the
/// variants carry.
///
/// `reader` is the outermost reader of the chain of the pass, lent and not
/// taken, so that whoever built the chain reads the counts of its filters
/// from it when this returns: those counts and
/// [`PerVarDistribs::num_vars`] are the `pass_stats` of a result in Python
/// and in TypeScript. The pass asks the reader for the genotypes alone.
///
/// Every statistic is calculated for each population over its individuals
/// alone, and a variant with too little data in a population is out of that
/// population's mean and in no bin of its histogram. A `config` that names
/// no statistic is a pass that reads the rows and gives a result of `None`
/// for each of the five, which the binding crates refuse before they get
/// here.
///
/// # Errors
///
/// A `poly_threshold` that is not a number from 0 to 1; what the reader
/// fails with; a block that holds no genotypes, which is a reader that was
/// asked for them and gave none, and a block of no variants, which is a
/// defect of a reader too; and a pass that gave no variant, whether its
/// source holds none or its steps kept none of them.
pub fn calc_per_var_distribs<R: BlockReader + ?Sized>(
    reader: &mut R,
    config: &PerVarDistribsConfig,
) -> Result<PerVarDistribs> {
    // A NaN is outside every range, so the comparison refuses it too.
    if !(0.0..=1.0).contains(&config.poly_threshold) {
        return Err(Error::PolyThresholdOutOfRange {
            value: config.poly_threshold,
        });
    }
    let asked = Asked::of(&config.stats);
    // The five statistics follow from the genotypes of a row, so no column
    // of a block is read and the reader is asked to fill none of them.
    reader.set_needs(Needs::GTS);
    let mut totals = Totals::of(config.pops.len(), config.bins.num_bins(), asked);
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        // The rows are cut out of the genotypes by the sizes the block
        // states, so those sizes are checked before anything is read.
        block.check()?;
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        if block.gts.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: Needs::GTS });
        }
        // `check` passed and the genotypes are not empty, so they are the
        // variants of the block times this number and it is one allele at
        // least: the rows are cut by it, and a cut of 0 is what the
        // standard library refuses with a panic.
        let alleles_per_var = block.alleles_per_var()?.max(1);
        add_the_block(&block, alleles_per_var, config, asked, &mut totals)?;
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so every one of them is a `u64`; and a pass of
        // more than 18446744073709551615 variants reads more rows than any
        // source holds.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
    }
    if num_vars == 0 {
        let filters = reader.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    Ok(the_distribs(&totals, config, asked, num_vars))
}

/// How many alleles one chunk of the pass holds: [`ROWS_PER_CHUNK`] rows of
/// `alleles_per_var` alleles, and one allele at least, because a cut of 0
/// is what the standard library refuses with a panic.
fn alleles_of_a_chunk(alleles_per_var: usize) -> usize {
    // A block of more alleles than a `usize` counts is refused before this,
    // and a chunk that saturated would be the whole block, which is a
    // chunking that gives the right numbers and no threads.
    ROWS_PER_CHUNK.saturating_mul(alleles_per_var).max(1)
}

/// It adds the statistics of every row of a block into `totals`.
///
/// Natively the chunks of rows are read on the threads of rayon, as section
/// 3 of `docs/architecture.md` asks: no row reads another, each chunk adds
/// up what its own rows give, and the chunks are added into `totals` in the
/// order of the block, so neither a count nor a sum depends on how many
/// threads there are. The threads are those of the pool the caller is
/// running in, and rayon's global pool only when the caller is in none.
///
/// # Errors
///
/// What the counts of one variant refuse: genotypes that are not a whole
/// number of genotypes of the ploidy, a variant of more alleles than a
/// count of them holds, an allele below the missing one, and an individual
/// of a population beyond the row. The error is the one of the first row
/// that has one, wherever the threads found it: which of two bad rows a
/// thread reaches first depends on how the chunks were shared out, and a
/// user who reports a damaged file has to get the same message every time,
/// so the rows are read again, one after another, to find the first.
#[cfg(not(target_family = "wasm"))]
fn add_the_block(
    block: &Block,
    alleles_per_var: usize,
    config: &PerVarDistribsConfig,
    asked: Asked,
    totals: &mut Totals,
) -> Result<()> {
    use rayon::iter::ParallelIterator;
    use rayon::slice::ParallelSlice;

    let num_pops = config.pops.len();
    let num_bins = config.bins.num_bins();
    let of_the_chunks: Result<Vec<Totals>> = block
        .gts
        .par_chunks(alleles_of_a_chunk(alleles_per_var))
        .map(|chunk| {
            let mut of_the_chunk = Totals::of(num_pops, num_bins, asked);
            add_the_rows(
                chunk,
                alleles_per_var,
                block.ploidy,
                config,
                asked,
                &mut of_the_chunk,
            )?;
            Ok(of_the_chunk)
        })
        .collect();
    match of_the_chunks {
        Ok(of_the_chunks) => {
            for of_the_chunk in &of_the_chunks {
                totals.add_the_chunk(of_the_chunk);
            }
            Ok(())
        }
        // The second pass costs a read of the block, and it is made only
        // where the block is refused and nothing of it is given.
        Err(of_a_thread) => {
            let mut read_again = Totals::of(num_pops, num_bins, asked);
            match add_the_chunks_one_by_one(block, alleles_per_var, config, asked, &mut read_again)
            {
                Err(of_the_first_row) => Err(of_the_first_row),
                // The rows are the same rows, so the second pass finds an
                // error too; the error of the threads is what is left if it
                // ever did not.
                Ok(()) => Err(of_a_thread),
            }
        }
    }
}

/// The same numbers, with the chunks read one after another, which is what
/// wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn add_the_block(
    block: &Block,
    alleles_per_var: usize,
    config: &PerVarDistribsConfig,
    asked: Asked,
    totals: &mut Totals,
) -> Result<()> {
    add_the_chunks_one_by_one(block, alleles_per_var, config, asked, totals)
}

/// The chunks of the block read one after another, each into accumulators
/// of its own that are added into `totals` before the next is read: what
/// wasm runs, and what the threads fall back on to find the first row that
/// is an error.
///
/// The chunks are the same chunks the threads read, and they are added in
/// the same order, so wasm and a native build add the values of a block up
/// in the same order and the addition of two floats is exact on every
/// machine popnei runs on.
///
/// # Errors
///
/// Those of [`add_the_block`], at the first row that has one.
fn add_the_chunks_one_by_one(
    block: &Block,
    alleles_per_var: usize,
    config: &PerVarDistribsConfig,
    asked: Asked,
    totals: &mut Totals,
) -> Result<()> {
    let mut of_the_chunk = Totals::of(config.pops.len(), config.bins.num_bins(), asked);
    for chunk in block.gts.chunks(alleles_of_a_chunk(alleles_per_var)) {
        of_the_chunk.forget_what_it_holds();
        add_the_rows(
            chunk,
            alleles_per_var,
            block.ploidy,
            config,
            asked,
            &mut of_the_chunk,
        )?;
        totals.add_the_chunk(&of_the_chunk);
    }
    Ok(())
}

/// It adds the statistics of every row of one chunk into `totals`, one
/// population after another.
///
/// `gts` holds whole rows of `alleles_per_var` alleles each. The counts of
/// one row over one population are taken once and the statistics that were
/// asked for follow from them: the observed heterozygosity from the counts
/// of the genotypes, the other four from the counts of the alleles.
///
/// # Errors
///
/// Those of [`add_the_block`], at the first row of the chunk that has one.
fn add_the_rows(
    gts: &[i8],
    alleles_per_var: usize,
    ploidy: usize,
    config: &PerVarDistribsConfig,
    asked: Asked,
    totals: &mut Totals,
) -> Result<()> {
    // One array of counts for every row and every population, which
    // `count_alleles_of` clears before it counts: a pass over a block
    // allocates nothing for a variant.
    let mut counts: AlleleCounts = [0; 128];
    for row in gts.chunks_exact(alleles_per_var) {
        for (pop, of_the_pop) in totals.pops.iter_mut().enumerate() {
            let individuals = config.pops.individuals(pop);
            // A population of every individual of the reader in its order
            // is counted by reading the row as it is, and the width of the
            // row says that it is the row of that reader: a `Pops` built
            // against another one would otherwise count individuals the
            // population does not hold.
            let of_the_whole_row = config.pops.is_all(pop)
                && individuals.len().saturating_mul(ploidy) == alleles_per_var;
            if asked.obs_het {
                let of_the_gts = if of_the_whole_row {
                    count_gts(row, ploidy)?
                } else {
                    count_gts_of(row, ploidy, individuals)?
                };
                if let Some(value) = config.obs_het.of_var(of_the_gts) {
                    of_the_pop.obs_het.add(value, &config.bins);
                }
            }
            if !asked.the_allele_counts() {
                continue;
            }
            let called_alleles = if of_the_whole_row {
                count_alleles(row, &mut counts)?
            } else {
                count_alleles_of(row, ploidy, individuals, &mut counts)?
            };
            if asked.the_maf()
                && let Some(value) = config.maf.of_var(&counts, called_alleles)
            {
                if asked.maf {
                    of_the_pop.maf.add(value, &config.bins);
                }
                if asked.poly_vars_ratio {
                    of_the_pop.poly.add(value, config.poly_threshold);
                }
            }
            if asked.exp_het
                && let Some(value) = config.exp_het.of_var(&counts, called_alleles, false)
            {
                of_the_pop.exp_het.add(value, &config.bins);
            }
            if asked.unbiased_exp_het
                && let Some(value) = config.exp_het.of_var(&counts, called_alleles, true)
            {
                of_the_pop.unbiased_exp_het.add(value, &config.bins);
            }
        }
    }
    Ok(())
}

/// The result of a pass, with the statistics that were asked for and
/// `None` for the ones that were not.
fn the_distribs(
    totals: &Totals,
    config: &PerVarDistribsConfig,
    asked: Asked,
    num_vars: u64,
) -> PerVarDistribs {
    let distrib = |asked_for: bool, of: fn(&OfAPop) -> &Accumulated| {
        asked_for.then(|| totals.distrib_of(&config.bins, of))
    };
    PerVarDistribs {
        obs_het: distrib(asked.obs_het, |pop| &pop.obs_het),
        maf: distrib(asked.maf, |pop| &pop.maf),
        exp_het: distrib(asked.exp_het, |pop| &pop.exp_het),
        unbiased_exp_het: distrib(asked.unbiased_exp_het, |pop| &pop.unbiased_exp_het),
        poly_vars_ratio: asked.poly_vars_ratio.then(|| totals.poly_vars_stats()),
        num_vars,
    }
}

#[cfg(test)]
mod pops {
    use super::{DEFAULT_POP_NAME, Pops};
    use crate::error::Error;

    /// The five diploid individuals of the worked example of
    /// `docs/specs/filters.md`, in the order of the source.
    fn the_five_individuals() -> Vec<String> {
        ["i1", "i2", "i3", "i4", "i5"]
            .iter()
            .map(|name| (*name).to_owned())
            .collect()
    }

    /// One population as `from_names` takes it, written as a user writes
    /// it.
    fn pop_of(name: &str, individuals: &[&str]) -> (String, Vec<String>) {
        (
            name.to_owned(),
            individuals
                .iter()
                .map(|individual| (*individual).to_owned())
                .collect(),
        )
    }

    /// The two populations of the worked example of "How it is verified"
    /// of the per variant distributions, pop1 of i1 and i2 and pop2 of i3,
    /// i4 and i5, which are the individuals 0 and 1 and the individuals 2,
    /// 3 and 4 of the source. Neither is every individual of the source.
    #[test]
    fn the_indices_of_the_two_populations_of_the_worked_example() {
        let pops = Pops::from_names(
            &[
                pop_of("pop1", &["i1", "i2"]),
                pop_of("pop2", &["i3", "i4", "i5"]),
            ],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.len(), 2);
        assert!(!pops.is_empty());
        assert_eq!(pops.name(0), "pop1");
        assert_eq!(pops.individuals(0), [0, 1]);
        assert!(!pops.is_all(0));
        assert_eq!(pops.name(1), "pop2");
        assert_eq!(pops.individuals(1), [2, 3, 4]);
        assert!(!pops.is_all(1));
    }

    /// The order of the populations is the one they were given in, and the
    /// order of the individuals of each the one the user named them in,
    /// which is not the order of the source here.
    #[test]
    fn the_populations_and_their_individuals_are_in_the_order_they_were_given() {
        let pops = Pops::from_names(
            &[pop_of("second", &["i5", "i1"]), pop_of("first", &["i3"])],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.name(0), "second");
        assert_eq!(pops.individuals(0), [4, 0]);
        assert_eq!(pops.name(1), "first");
        assert_eq!(pops.individuals(1), [2]);
    }

    /// The name a user wrote and the population it is in are what they
    /// have to look at, so the error carries both.
    #[test]
    fn a_name_that_is_not_an_individual_is_refused_with_its_population() {
        let error = Pops::from_names(&[pop_of("pop1", &["i1", "nope"])], &the_five_individuals())
            .unwrap_err();
        assert!(
            matches!(&error, Error::IndividualOfAPopNotInThePass { pop, name }
                if pop == "pop1" && name == "nope"),
            "{error:?}"
        );
    }

    /// pyNei counts an individual named twice twice, and popnei refuses it:
    /// a count that is wrong and says nothing is what popnei never gives.
    #[test]
    fn a_name_twice_in_one_population_is_refused_with_its_population() {
        let error = Pops::from_names(
            &[pop_of("pop1", &["i1", "i2", "i1"])],
            &the_five_individuals(),
        )
        .unwrap_err();
        assert!(
            matches!(&error, Error::IndividualNamedTwiceInAPop { pop, name }
                if pop == "pop1" && name == "i1"),
            "{error:?}"
        );
    }

    /// A population with no individual has no value for any statistic,
    /// where pyNei gives NaN, so the name of that population is the error.
    #[test]
    fn a_population_with_no_individual_is_refused_with_its_name() {
        let error = Pops::from_names(
            &[pop_of("pop1", &["i1"]), pop_of("empty", &[])],
            &the_five_individuals(),
        )
        .unwrap_err();
        assert!(
            matches!(&error, Error::PopWithNoIndividual { pop } if pop == "empty"),
            "{error:?}"
        );
    }

    /// `pops` with no population at all would leave a result with nothing
    /// in it, where pyNei gives one with no column.
    #[test]
    fn pops_with_no_population_is_refused() {
        let error = Pops::from_names(&[], &the_five_individuals()).unwrap_err();
        assert!(matches!(&error, Error::NoPop), "{error:?}");
    }

    /// An individual that is in two populations is taken, as in pyNei,
    /// whose `test_maf_stats` names one individual in both of its
    /// populations. Each population holds it once.
    #[test]
    fn an_individual_in_two_populations_is_taken() {
        let pops = Pops::from_names(
            &[pop_of("pop1", &["i1", "i3"]), pop_of("pop2", &["i3", "i4"])],
            &the_five_individuals(),
        )
        .unwrap();
        assert_eq!(pops.individuals(0), [0, 2]);
        assert_eq!(pops.individuals(1), [2, 3]);
    }

    /// What a calculation that was given no `pops` works over: one
    /// population of every individual of the reader, in its order, named
    /// as pyNei names it.
    #[test]
    fn all_gives_one_population_named_pop_of_every_individual() {
        let pops = Pops::all(5);
        assert_eq!(pops.len(), 1);
        assert!(!pops.is_empty());
        assert_eq!(pops.name(0), DEFAULT_POP_NAME);
        assert_eq!(pops.name(0), "pop");
        assert_eq!(pops.individuals(0), [0, 1, 2, 3, 4]);
        assert!(pops.is_all(0));
    }

    /// A population a user named that turns out to be every individual in
    /// the order of the source is every individual: what `is_all` tells a
    /// caller is what the indices are, and not which constructor made
    /// them. One in another order, or one that leaves an individual out,
    /// is not.
    #[test]
    fn a_population_of_every_individual_in_the_order_of_the_source_is_all_of_them() {
        let pops = Pops::from_names(
            &[
                pop_of("in order", &["i1", "i2", "i3", "i4", "i5"]),
                pop_of("another order", &["i1", "i2", "i3", "i5", "i4"]),
                pop_of("four of them", &["i1", "i2", "i3", "i4"]),
            ],
            &the_five_individuals(),
        )
        .unwrap();
        assert!(pops.is_all(0));
        assert!(!pops.is_all(1));
        assert!(!pops.is_all(2));
    }
}

#[cfg(test)]
mod fixtures {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::{Path, PathBuf};

    use crate::block::{Block, BlockReader};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{
        AlleleCounts, GtCounts, Needs, count_alleles_of, count_gts, count_gts_of,
    };

    /// The individuals of pop1 of the worked example of "How it is
    /// verified" of the per variant distributions, i1 and i2, as indices
    /// among the five individuals of the source.
    pub(super) const POP1: [usize; 2] = [0, 1];

    /// The individuals of pop2 of that worked example, i3, i4 and i5.
    pub(super) const POP2: [usize; 3] = [2, 3, 4];

    /// The six variants of five diploid individuals of the worked example
    /// of the per variant distributions, which is the one of
    /// `docs/specs/filters.md`, one row of genotypes each:
    /// `0/0 0/1 0/0 0/0 0/.`, `0/0 0/1 0/0 ./. 0/.`, `0/1 2/3 0/1 2/3 ./.`,
    /// five missing genotypes, `0/0 0/0 0/0 0/0 1/1` and
    /// `0/. ./. ./. ./. ./.`.
    pub(super) const THE_SIX_VARIANTS: [[i8; 10]; 6] = [
        [0, 0, 0, 1, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 1, 0, 0, -1, -1, 0, -1],
        [0, 1, 2, 3, 0, 1, 2, 3, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ];

    /// The three variants of five diploid individuals of the worked example
    /// of the expected heterozygosity: `0/0 2/1 0/0 0/0 0/.`,
    /// `0/0 0/0 0/1 1/0 ./.` and five missing genotypes.
    pub(super) const THE_THREE_VARIANTS: [[i8; 10]; 3] = [
        [0, 0, 2, 1, 0, 0, 0, 0, 0, -1],
        [0, 0, 0, 0, 0, 1, 1, 0, -1, -1],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ];

    /// The counts of the genotypes of one variant of a worked example over
    /// the individuals of one population, which the statistic then reads.
    pub(super) fn gt_counts_of(gts: &[i8], pop: &[usize]) -> GtCounts {
        count_gts_of(gts, 2, pop).expect("the genotype counts of the population")
    }

    /// The counts of the alleles of one variant of a worked example over
    /// the individuals of one population, with their sum, the called
    /// alleles of the population.
    pub(super) fn allele_counts_of(gts: &[i8], pop: &[usize]) -> (AlleleCounts, u32) {
        let mut counts: AlleleCounts = [0; 128];
        let called_alleles = count_alleles_of(gts, 2, pop, &mut counts)
            .expect("the allele counts of the population");
        (counts, called_alleles)
    }

    /// An array of allele counts with `counts[a]` of the allele a, which is
    /// how a variant of the panel and of `many.vcf` is written into a test
    /// from the numbers the reference program printed.
    pub(super) fn allele_counts(of_each_allele: &[u32]) -> AlleleCounts {
        let mut counts: AlleleCounts = [0; 128];
        for (entry, count) in counts.iter_mut().zip(of_each_allele) {
            *entry = *count;
        }
        counts
    }

    /// The individuals of `popA` of `many.vcf`, the first 20 of the 50,
    /// which `tests/reference/stats/many_pops.txt` names.
    pub(super) fn pop_a() -> Vec<usize> {
        (0..20).collect()
    }

    /// The individuals of `popB` of `many.vcf`, the other 30.
    pub(super) fn pop_b() -> Vec<usize> {
        (20..50).collect()
    }

    /// The reference VCFs live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    fn many_vcf() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/vcf/many.vcf")
    }

    /// The first three variants of `many.vcf`, the 500 variants of 50
    /// diploid individuals of `docs/specs/io_vcf.md`, at the positions
    /// 1000, 1037 and 1074, with every variant given.
    pub(super) fn the_first_block_of_many_vcf() -> Block {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(3),
        };
        let mut reader =
            VcfReader::<BufReader<File>>::from_path(&many_vcf(), options).expect("many.vcf");
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let block = reader
            .next_block()
            .expect("the first block of many.vcf")
            .expect("many.vcf has variants");
        assert_eq!(block.num_vars, 3);
        block
    }

    /// The genotypes of one variant of a block, which a statistic counts.
    pub(super) fn gts_of(block: &Block, var: usize) -> &[i8] {
        block.variant(var).expect("the variant of the block").gts()
    }

    /// The counts of the genotypes of one variant over every individual of
    /// the row, in the order of the source.
    pub(super) fn gt_counts(gts: &[i8]) -> GtCounts {
        count_gts(gts, 2).expect("the genotype counts of the variant")
    }

    /// That a statistic of one variant has the value the spec gives, within
    /// `tolerance`, which is one unit of the last digit the reference
    /// program or the spec prints.
    pub(super) fn assert_value(found: Option<f64>, expected: f64, tolerance: f64, what: &str) {
        let Some(found) = found else {
            panic!("{what} has no value, and it is {expected}");
        };
        assert!(
            (found - expected).abs() <= tolerance,
            "{what} is {found}, and it is {expected} within {tolerance}"
        );
    }
}

#[cfg(test)]
mod hist {
    use super::{HistBins, MAX_NUM_BINS};
    use crate::error::Error;

    /// One unit of the last digit that numpy and pyNei print of a
    /// logarithmic edge, which is the relative difference the spec allows
    /// between the two libraries there: `powf` and `ln` are not rounded the
    /// same on every platform, and popnei promises no more than this for
    /// them.
    const OF_A_LOGARITHMIC_EDGE: f64 = 1e-12;

    /// The bins of the worked example of the per variant distributions,
    /// four of equal width from 0 to 1, whose edges the spec gives. They
    /// are compared exactly: `start + i * step` uses the four operations
    /// alone, which are rounded the same everywhere.
    #[test]
    fn the_edges_of_four_linear_bins_from_0_to_1() {
        let bins = HistBins::linear(0.0, 1.0, 4).unwrap();
        assert_eq!(bins.num_bins(), 4);
        assert_eq!(bins.edges(), [0.0, 0.25, 0.5, 0.75, 1.0]);
    }

    /// The default histogram of the per variant distributions, 40 bins of
    /// equal width from 0 to 1. The literals are what `numpy.linspace(0, 1,
    /// 41)` gives, printed by Python at the shortest text that reads back
    /// as the same float64, and they are compared exactly: a value that
    /// falls on an edge has to fall on the same side of it in popnei and in
    /// pyNei, or a count of the histogram is one out.
    #[test]
    fn the_edges_of_the_forty_default_linear_bins() {
        let bins = HistBins::linear(0.0, 1.0, 40).unwrap();
        assert_eq!(bins.num_bins(), 40);
        assert_eq!(
            bins.edges(),
            [
                0.0,
                0.025,
                0.05,
                0.07500000000000001,
                0.1,
                0.125,
                0.15000000000000002,
                0.17500000000000002,
                0.2,
                0.225,
                0.25,
                0.275,
                0.30000000000000004,
                0.325,
                0.35000000000000003,
                0.375,
                0.4,
                0.42500000000000004,
                0.45,
                0.47500000000000003,
                0.5,
                0.525,
                0.55,
                0.5750000000000001,
                0.6000000000000001,
                0.625,
                0.65,
                0.675,
                0.7000000000000001,
                0.7250000000000001,
                0.75,
                0.775,
                0.8,
                0.8250000000000001,
                0.8500000000000001,
                0.875,
                0.9,
                0.925,
                0.9500000000000001,
                0.9750000000000001,
                1.0,
            ]
        );
    }

    /// The example of `test_logarithmic_bins_span_the_given_range` of
    /// pyNei's `test_hist.py`: four bins of equal ratio from 0.01 to 100
    /// have the edges 0.01, 0.1, 1, 10 and 100. They go through `log10` and
    /// a power of 10, which are not rounded the same on every platform, so
    /// they are compared within one unit of the last digit printed.
    #[test]
    fn the_edges_of_the_logarithmic_example_of_pynei() {
        let bins = HistBins::logarithmic(0.01, 100.0, 4).unwrap();
        assert_eq!(bins.num_bins(), 4);
        for (found, expected) in bins.edges().iter().zip([0.01, 0.1, 1.0, 10.0, 100.0]) {
            assert!(
                (found - expected).abs() <= OF_A_LOGARITHMIC_EDGE * expected,
                "the edge is {found}, and it is {expected}"
            );
        }
    }

    /// A value on an edge inside the range is in the bin that starts there
    /// and not in the one that ends there, which is what `numpy.histogram`
    /// does: every bin but the last takes its left edge and not its right
    /// one.
    #[test]
    fn a_value_on_an_interior_edge_is_in_the_bin_that_starts_there() {
        let bins = HistBins::linear(0.0, 1.0, 4).unwrap();
        assert_eq!(bins.bin_of(0.25), Some(1));
        assert_eq!(bins.bin_of(0.5), Some(2));
        assert_eq!(bins.bin_of(0.75), Some(3));
        // And a value between two edges is in the bin they are the ends of.
        assert_eq!(bins.bin_of(0.0), Some(0));
        assert_eq!(bins.bin_of(0.1), Some(0));
        assert_eq!(bins.bin_of(0.3), Some(1));
        assert_eq!(bins.bin_of(0.99), Some(3));
    }

    /// The last bin takes its right edge too, so an observed heterozygosity
    /// or a maf of exactly 1 is counted, in the last bin. The bins of equal
    /// ratio do the same at their own end.
    #[test]
    fn a_value_on_the_last_edge_is_in_the_last_bin() {
        assert_eq!(HistBins::linear(0.0, 1.0, 4).unwrap().bin_of(1.0), Some(3));
        assert_eq!(
            HistBins::linear(0.0, 1.0, 40).unwrap().bin_of(1.0),
            Some(39)
        );
        let of_equal_ratio = HistBins::logarithmic(0.01, 100.0, 4).unwrap();
        assert_eq!(of_equal_ratio.bin_of(100.0), Some(3));
        assert_eq!(of_equal_ratio.bin_of(0.01), Some(0));
    }

    /// A value outside the range is in no bin, and it still counts in the
    /// mean: pyNei's unbiased expected heterozygosity above 1 is the case
    /// the spec names.
    #[test]
    fn a_value_outside_the_range_is_in_no_bin() {
        let bins = HistBins::linear(0.0, 1.0, 4).unwrap();
        assert_eq!(bins.bin_of(-0.1), None);
        assert_eq!(bins.bin_of(1.1), None);
        assert_eq!(bins.bin_of(1.1685), None);
        assert_eq!(bins.bin_of(f64::NAN), None);
        assert_eq!(bins.bin_of(f64::INFINITY), None);
        assert_eq!(bins.bin_of(f64::NEG_INFINITY), None);
        let of_equal_ratio = HistBins::logarithmic(0.01, 100.0, 4).unwrap();
        assert_eq!(of_equal_ratio.bin_of(0.0), None);
        assert_eq!(of_equal_ratio.bin_of(1000.0), None);
    }

    /// A histogram of no bin counts nothing, so it is refused instead of
    /// built: numpy gives one edge for it and a histogram no value falls
    /// in.
    #[test]
    fn a_histogram_of_no_bin_is_refused() {
        for built in [
            HistBins::linear(0.0, 1.0, 0),
            HistBins::logarithmic(0.01, 1.0, 0),
        ] {
            let error = built.unwrap_err();
            assert!(matches!(&error, Error::HistWithNoBin), "{error:?}");
        }
    }

    /// The range has to run from a number up to a larger one: two ends the
    /// wrong way round, two that are equal, and one that is NaN or
    /// infinite, which would leave every edge between them NaN.
    #[test]
    fn a_range_that_does_not_go_up_is_refused() {
        for (start, end) in [
            (1.0, 0.0),
            (0.5, 0.5),
            (f64::NAN, 1.0),
            (0.0, f64::NAN),
            (f64::NEG_INFINITY, f64::INFINITY),
            (0.0, f64::INFINITY),
        ] {
            let error = HistBins::linear(start, end, 4).unwrap_err();
            assert!(
                matches!(&error, Error::HistRangeNotGoingUp { .. }),
                "{start} to {end}: {error:?}"
            );
            let error = HistBins::logarithmic(start, end, 4).unwrap_err();
            assert!(
                matches!(&error, Error::HistRangeNotGoingUp { .. }),
                "{start} to {end}: {error:?}"
            );
        }
    }

    /// The distance between the two ends of the range is a number too: the
    /// width of a bin is that distance over the bins, and 4 bins from
    /// -1e308 to 1e308 have the width infinity and the edges NaN,
    /// infinity, infinity, infinity and 1e308, which do not go up, so the
    /// search for the bin of a value puts every value in the first bin.
    /// numpy refuses the same range with "Too many bins for data range.
    /// Cannot create 4 finite-sized bins."
    #[test]
    fn a_range_whose_two_ends_are_further_apart_than_a_float64_goes_is_refused() {
        for (start, end) in [(-1e308, 1e308), (f64::MIN, f64::MAX)] {
            let error = HistBins::linear(start, end, 4).unwrap_err();
            assert!(
                matches!(&error, Error::HistRangeTooWide { .. }),
                "{start} to {end}: {error:?}"
            );
        }
        // No range of bins of equal ratio has an infinite width: its two
        // ends are above 0, so the distance between them is the end at
        // most, and the widths are taken over the base 10 logarithms of the
        // ends, which lie between -324 and 309.
        let of_equal_ratio =
            HistBins::logarithmic(1e-300, 1e300, 4).expect("the widest range of equal ratios");
        assert!(of_equal_ratio.edges().iter().all(|edge| edge.is_finite()));
    }

    /// A histogram has [`MAX_NUM_BINS`] bins at most, whatever a user
    /// writes for `num_bins`: the counts of one of 2^60 bins are a vector
    /// no machine gives, which is a panic of the allocation and, in Python,
    /// a `PanicException` that derives from `BaseException` and ends a
    /// notebook.
    #[test]
    fn a_histogram_of_more_bins_than_a_person_reads_is_refused() {
        for num_bins in [MAX_NUM_BINS.saturating_add(1), 1 << 40, usize::MAX] {
            for built in [
                HistBins::linear(0.0, 1.0, num_bins),
                HistBins::logarithmic(0.01, 100.0, num_bins),
            ] {
                let error = built.unwrap_err();
                assert!(
                    matches!(&error, Error::HistTooManyBins { num_bins: found, largest }
                        if *found == num_bins && *largest == MAX_NUM_BINS),
                    "{num_bins}: {error:?}"
                );
            }
        }
        // The bound itself is built, and its bins are counted.
        assert_eq!(
            HistBins::linear(0.0, 1.0, MAX_NUM_BINS)
                .expect("the largest histogram")
                .num_bins(),
            MAX_NUM_BINS
        );
    }

    /// Bins of equal ratio start above 0, since each edge is the one before
    /// it times a fixed factor and no factor takes 0 anywhere. The two
    /// ranges are those of `test_logarithmic_bins_refuse_a_non_positive_range`
    /// of pyNei's `test_hist.py`, which refuses both as well.
    #[test]
    fn a_logarithmic_range_that_starts_at_0_or_below_is_refused() {
        for start in [0.0, -1.0] {
            let error = HistBins::logarithmic(start, 10.0, 4).unwrap_err();
            assert!(
                matches!(&error, Error::HistLogRangeNotAboveZero { .. }),
                "{start}: {error:?}"
            );
            // The same range is a range of bins of equal width.
            assert!(HistBins::linear(start, 10.0, 4).is_ok());
        }
    }
}

#[cfg(test)]
mod obs_het {
    use super::ObsHet;
    use super::fixtures::{
        POP1, POP2, THE_SIX_VARIANTS, assert_value, gt_counts, gt_counts_of, gts_of,
        the_first_block_of_many_vcf,
    };
    use crate::variant::GtCounts;

    /// One unit of the last of the six digits plink2 prints of an observed
    /// heterozygosity, which is the tolerance the spec gives for every
    /// comparison of this statistic with a reference program.
    const OF_A_PRINTED_VALUE: f64 = 1e-6;

    /// The observed heterozygosity of the six variants of the worked
    /// example of the per variant distributions, in pop1, i1 and i2, and in
    /// pop2, i3, i4 and i5, with `min_num_individuals` 1. The values are
    /// those of `_calc_obs_het_per_var` of pyNei at commit ef0ca6e on these
    /// genotypes, as the spec's table gives them.
    ///
    /// Variant 1 in pop2 has one called allele from its half called
    /// genotype, which this statistic does not count: 2 called genotypes,
    /// none heterozygous. Variant 6 in pop1 has no called genotype at all.
    #[test]
    fn the_values_of_the_worked_example_of_the_pass() {
        let obs_het = ObsHet::new(1);
        let of_pop1 = [Some(0.5), Some(0.5), Some(1.0), None, Some(0.0), None];
        let of_pop2 = [Some(0.0), Some(0.0), Some(1.0), None, Some(0.0), None];
        for (var, gts) in THE_SIX_VARIANTS.iter().enumerate() {
            for (pop, expected) in [(&POP1[..], of_pop1[var]), (&POP2[..], of_pop2[var])] {
                let found = obs_het.of_var(gt_counts_of(gts, pop));
                match expected {
                    Some(expected) => assert_value(
                        found,
                        expected,
                        OF_A_PRINTED_VALUE,
                        &format!("the observed heterozygosity of variant {var}"),
                    ),
                    None => assert_eq!(found, None, "variant {var} has a value"),
                }
            }
        }
    }

    /// The literals of `var0000` of the panel, which the `--hardy` report
    /// of plink2 v2.0.0-a.7.7 prints per population: the genotypes
    /// homozygous for the reference allele, heterozygous and homozygous for
    /// the other, and `O(HET_A1)`, the heterozygous over the three, which
    /// for a biallelic variant is this statistic.
    #[test]
    fn the_literals_of_var0000_of_the_panel() {
        let obs_het = ObsHet::new(20);
        // p0 has 48 individuals and all of them are called here, p1 68 with
        // one missing, p2 84 with one missing, and the panel 200 with two.
        for (called, missing, het, expected, pop) in [
            (48_u32, 0_u32, 8_u32, 0.166_667, "p0"),
            (67, 1, 34, 0.507_463, "p1"),
            (83, 1, 2, 0.024_096_4, "p2"),
            (198, 2, 44, 0.222_222, "the 200 individuals"),
        ] {
            let counts = GtCounts {
                called,
                missing,
                het,
            };
            assert_value(
                obs_het.of_var(counts),
                expected,
                OF_A_PRINTED_VALUE,
                &format!("the observed heterozygosity of var0000 in {pop}"),
            );
        }
    }

    /// The first three variants of `many.vcf`, over its 50 individuals,
    /// which have a half called genotype and a third allele among them.
    /// bcftools 1.24 counts the heterozygous genotypes and the called ones
    /// of each, which `tests/reference/stats/many.counts.tsv` holds: 9 of
    /// 46, 24 of 47 and 32 of 47 at the positions 1000, 1037 and 1074. A
    /// heterozygous genotype is one whose alleles differ, whichever they
    /// are, and a half called one is missing.
    #[test]
    fn the_first_three_variants_of_many_vcf() {
        let obs_het = ObsHet::new(5);
        let block = the_first_block_of_many_vcf();
        let expected = [
            (1000_u64, 9_u32, 46_u32, 0.195_652),
            (1037, 24, 47, 0.510_638),
            (1074, 32, 47, 0.680_851),
        ];
        for (var, (pos, het, called, value)) in expected.into_iter().enumerate() {
            let variant = block.variant(var).expect("the variant");
            assert_eq!(variant.pos(), Some(pos));
            let counts = gt_counts(gts_of(&block, var));
            assert_eq!((counts.het, counts.called), (het, called));
            assert_value(
                obs_het.of_var(counts),
                value,
                OF_A_PRINTED_VALUE,
                &format!("the observed heterozygosity at {pos}"),
            );
        }
    }

    /// A population with fewer called genotypes than `min_num_individuals`
    /// has no value, and one with exactly that many keeps it: pyNei exempts
    /// this statistic from the threshold and popnei holds it to it. The
    /// test is on called genotypes, whole ones, where the other statistics
    /// count called alleles over the ploidy.
    #[test]
    fn a_population_below_min_num_individuals_has_no_value() {
        let obs_het = ObsHet::new(20);
        let of = |called| {
            obs_het.of_var(GtCounts {
                called,
                missing: 1,
                het: 5,
            })
        };
        assert_eq!(of(19), None);
        assert_value(of(20), 0.25, OF_A_PRINTED_VALUE, "20 called genotypes");
        // A threshold of 0 keeps every population that called something.
        assert_value(
            ObsHet::new(0).of_var(GtCounts {
                called: 1,
                missing: 4,
                het: 1,
            }),
            1.0,
            OF_A_PRINTED_VALUE,
            "one called genotype at a threshold of 0",
        );
    }

    /// A population that has called nothing at a variant has no value,
    /// whatever the threshold: the share of its called genotypes that are
    /// heterozygous is 0 over 0.
    #[test]
    fn a_population_with_nothing_called_has_no_value() {
        let nothing_called = GtCounts {
            called: 0,
            missing: 5,
            het: 0,
        };
        assert_eq!(ObsHet::new(0).of_var(nothing_called), None);
        assert_eq!(ObsHet::new(20).of_var(nothing_called), None);
    }
}

#[cfg(test)]
mod maf {
    use super::Maf;
    use super::fixtures::{
        POP1, POP2, THE_SIX_VARIANTS, allele_counts, allele_counts_of, assert_value, gts_of, pop_a,
        pop_b, the_first_block_of_many_vcf,
    };
    use crate::error::Error;
    use crate::variant::count_alleles_of;

    /// One unit of the last of the six digits plink2 prints of a frequency,
    /// which is the tolerance the spec gives for every comparison of this
    /// statistic with a reference program.
    const OF_A_PRINTED_VALUE: f64 = 1e-6;

    /// The major allele frequency of the six variants of the worked example
    /// of the per variant distributions, in pop1, i1 and i2, and in pop2,
    /// i3, i4 and i5, with `min_num_individuals` 1. The values are those of
    /// `_calc_maf_per_var` of pyNei at commit ef0ca6e on these genotypes,
    /// as the spec's table gives them.
    ///
    /// Variant 1 in pop2 has 5 called alleles, the fifth from its half
    /// called genotype, and all of them are the allele 0. Variant 6 in pop1
    /// has one called allele, half a genotype, below the threshold of 1.
    #[test]
    fn the_values_of_the_worked_example_of_the_pass() {
        let maf = Maf::new(2, 1).unwrap();
        let of_pop1 = [Some(0.75), Some(0.75), Some(0.25), None, Some(1.0), None];
        let of_pop2 = [
            Some(1.0),
            Some(1.0),
            Some(0.25),
            None,
            Some(0.666_667),
            None,
        ];
        for (var, gts) in THE_SIX_VARIANTS.iter().enumerate() {
            for (pop, expected) in [(&POP1[..], of_pop1[var]), (&POP2[..], of_pop2[var])] {
                let (counts, called_alleles) = allele_counts_of(gts, pop);
                let found = maf.of_var(&counts, called_alleles);
                match expected {
                    Some(expected) => assert_value(
                        found,
                        expected,
                        OF_A_PRINTED_VALUE,
                        &format!("the maf of variant {var}"),
                    ),
                    None => assert_eq!(found, None, "variant {var} has a value"),
                }
            }
        }
    }

    /// The literals of `var0000` of the panel, from the `--freq` report of
    /// plink2 v2.0.0-a.7.7, which prints the called alleles of each
    /// population and the frequency of each of the two alleles. The counts
    /// of the two come from the genotypes of the `--hardy` report of the
    /// same run: in p0, 39 homozygotes for the reference allele, 8
    /// heterozygotes and 1 homozygote for the other are 39 * 2 + 8 = 86 of
    /// the reference allele and 8 + 1 * 2 = 10 of the other. In p1 the
    /// alternative allele is the major one.
    #[test]
    fn the_literals_of_var0000_of_the_panel() {
        let maf = Maf::new(2, 20).unwrap();
        for (of_each_allele, called_alleles, expected, pop) in [
            ([86_u32, 10_u32], 96_u32, 0.895_833, "p0"),
            ([64, 70], 134, 0.522_388, "p1"),
            ([164, 2], 166, 0.987_952, "p2"),
            ([314, 82], 396, 0.792_929, "the 200 individuals"),
        ] {
            assert_eq!(of_each_allele[0] + of_each_allele[1], called_alleles);
            assert_value(
                maf.of_var(&allele_counts(&of_each_allele), called_alleles),
                expected,
                OF_A_PRINTED_VALUE,
                &format!("the maf of var0000 in {pop}"),
            );
        }
    }

    /// The first three variants of `many.vcf`, in `popA`, its first 20
    /// individuals, and in `popB`, the other 30, which have a half called
    /// genotype and a third allele among them. bcftools 1.24 counts the
    /// called alleles of each population and each alternative allele, which
    /// `tests/reference/stats/many.counts.tsv` holds, and the spec reads
    /// the largest count of each: 32 of 36, 27 of 39 and 21 of 40 in
    /// `popA`, and 52 of 57, 35 of 56 and 24 of 55 in `popB`. A half called
    /// genotype gives its called allele.
    #[test]
    fn the_first_three_variants_of_many_vcf() {
        let maf = Maf::new(2, 0).unwrap();
        let block = the_first_block_of_many_vcf();
        let expected = [
            (1000_u64, [36_u32, 57_u32], [0.888_889, 0.912_281]),
            (1037, [39, 56], [0.692_308, 0.625]),
            (1074, [40, 55], [0.525, 0.436_364]),
        ];
        for (var, (pos, of_each_pop, values)) in expected.into_iter().enumerate() {
            let variant = block.variant(var).expect("the variant");
            assert_eq!(variant.pos(), Some(pos));
            for (pop, (individuals, name)) in [(pop_a(), "popA"), (pop_b(), "popB")]
                .into_iter()
                .enumerate()
            {
                let mut counts = [0_u32; 128];
                let called_alleles =
                    count_alleles_of(gts_of(&block, var), 2, &individuals, &mut counts).unwrap();
                assert_eq!(called_alleles, of_each_pop[pop]);
                assert_value(
                    maf.of_var(&counts, called_alleles),
                    values[pop],
                    OF_A_PRINTED_VALUE,
                    &format!("the maf at {pos} in {name}"),
                );
            }
        }
    }

    /// A population that has called fewer than `min_num_individuals`
    /// genotypes has no value, and one with exactly that many keeps it. The
    /// test is on the called alleles against that many times the ploidy,
    /// which is the same test pyNei makes by dividing, and it is the ploidy
    /// the statistic was built with and not the count of the alleles.
    #[test]
    fn a_population_below_min_num_individuals_has_no_value() {
        let of_diploids = Maf::new(2, 20).unwrap();
        let counts = allele_counts(&[30, 10]);
        assert_eq!(of_diploids.of_var(&counts, 39), None);
        assert_value(
            of_diploids.of_var(&counts, 40),
            0.75,
            OF_A_PRINTED_VALUE,
            "20 called genotypes of a diploid population",
        );
        // The same 40 called alleles are 10 genotypes of a tetraploid
        // population, below a threshold of 20 and above one of 10.
        let of_tetraploids = Maf::new(4, 20).unwrap();
        assert_eq!(of_tetraploids.of_var(&counts, 40), None);
        assert_value(
            Maf::new(4, 10).unwrap().of_var(&counts, 40),
            0.75,
            OF_A_PRINTED_VALUE,
            "10 called genotypes of a tetraploid population",
        );
    }

    /// A population that has called nothing at a variant has no value,
    /// whatever the threshold: its commonest allele was called 0 times out
    /// of 0.
    #[test]
    fn a_population_with_nothing_called_has_no_value() {
        let nothing_called = allele_counts(&[]);
        assert_eq!(Maf::new(2, 0).unwrap().of_var(&nothing_called, 0), None);
        assert_eq!(Maf::new(2, 20).unwrap().of_var(&nothing_called, 0), None);
    }

    /// A ploidy of 0 would ask a population for 0 called alleles and turn
    /// the `min_num_individuals` test off, and one above 255 is above the
    /// largest ploidy a reader of popnei gives, so the constructor refuses
    /// both instead of building a statistic whose values say nothing.
    #[test]
    fn the_constructor_refuses_a_ploidy_of_0_and_of_256() {
        for ploidy in [0, 256, usize::MAX] {
            let error = Maf::new(ploidy, 20).unwrap_err();
            assert!(
                matches!(&error, Error::StatPloidyOutOfRange { kind, value, largest }
                    if *kind == "ploidy" && *value == ploidy && *largest == 255),
                "{ploidy}: {error:?}"
            );
        }
        assert!(Maf::new(1, 20).is_ok());
        assert!(Maf::new(255, 20).is_ok());
    }
}

#[cfg(test)]
mod exp_het {
    use super::ExpHet;
    use super::fixtures::{
        POP1, POP2, THE_SIX_VARIANTS, THE_THREE_VARIANTS, allele_counts, allele_counts_of,
        assert_value,
    };
    use crate::error::Error;

    /// One unit of the last of the six digits plink2 and the spec print of
    /// an expected heterozygosity, which is the tolerance the spec gives
    /// for every comparison of this statistic with a reference program.
    const OF_A_PRINTED_VALUE: f64 = 1e-6;

    /// The unbiased expected heterozygosity of the six variants of the
    /// worked example of the per variant distributions, in pop1, i1 and i2,
    /// and in pop2, i3, i4 and i5, with `min_num_individuals` 1. The values
    /// are those of `_calc_unbiased_exp_het_per_var` of pyNei at commit
    /// ef0ca6e on these genotypes, as the spec's table gives them; at the
    /// ploidy 2 of the data popnei's correction and pyNei's are the same
    /// factor.
    #[test]
    fn the_unbiased_values_of_the_worked_example_of_the_pass() {
        let exp_het = ExpHet::new(2, 2, 1).unwrap();
        let of_pop1 = [Some(0.5), Some(0.5), Some(1.0), None, Some(0.0), None];
        let of_pop2 = [Some(0.0), Some(0.0), Some(1.0), None, Some(0.533_333), None];
        for (var, gts) in THE_SIX_VARIANTS.iter().enumerate() {
            for (pop, expected) in [(&POP1[..], of_pop1[var]), (&POP2[..], of_pop2[var])] {
                let (counts, called_alleles) = allele_counts_of(gts, pop);
                let found = exp_het.of_var(&counts, called_alleles, true);
                match expected {
                    Some(expected) => assert_value(
                        found,
                        expected,
                        OF_A_PRINTED_VALUE,
                        &format!("the unbiased expected heterozygosity of variant {var}"),
                    ),
                    None => assert_eq!(found, None, "variant {var} has a value"),
                }
            }
        }
    }

    /// The three variants of the worked example of the expected
    /// heterozygosity, in pop1, i1 and i2, and in pop2, i3, i4 and i5, with
    /// `min_num_individuals` 1, plain and unbiased. They are the numbers
    /// `test_calc_exp_het` of pyNei's `test/test_diversity.py` asserts.
    ///
    /// Variant 1 in pop2 has 5 called alleles, four from two genotypes and
    /// one from the half called one, all of the allele 0, so both values
    /// are 0.
    #[test]
    fn the_values_of_the_worked_example_of_the_expected_heterozygosity() {
        let exp_het = ExpHet::new(2, 2, 1).unwrap();
        let plain = [
            [Some(0.625), Some(0.0)],
            [Some(0.0), Some(0.5)],
            [None, None],
        ];
        let unbiased = [
            [Some(0.833_333), Some(0.0)],
            [Some(0.0), Some(0.666_667)],
            [None, None],
        ];
        for (var, gts) in THE_THREE_VARIANTS.iter().enumerate() {
            for (pop, individuals) in [(0, &POP1[..]), (1, &POP2[..])] {
                let (counts, called_alleles) = allele_counts_of(gts, individuals);
                for (unbiased_one, expected) in
                    [(false, plain[var][pop]), (true, unbiased[var][pop])]
                {
                    let found = exp_het.of_var(&counts, called_alleles, unbiased_one);
                    match expected {
                        Some(expected) => assert_value(
                            found,
                            expected,
                            OF_A_PRINTED_VALUE,
                            &format!("variant {var} of pop {pop}, unbiased {unbiased_one}"),
                        ),
                        None => assert_eq!(found, None, "variant {var} of pop {pop} has a value"),
                    }
                }
            }
        }
    }

    /// The plain literals of `var0000` of the panel, which the `--hardy`
    /// report of plink2 v2.0.0-a.7.7 prints as `E(HET_A1)`, the frequency of
    /// heterozygotes expected under Hardy Weinberg, which for a biallelic
    /// variant is 1 - p² - q². The counts of the two alleles come from the
    /// genotypes of the same report, as the maf's test of the panel says.
    #[test]
    fn the_plain_literals_of_var0000_of_the_panel() {
        let exp_het = ExpHet::new(2, 2, 20).unwrap();
        for (of_each_allele, called_alleles, expected, pop) in [
            ([86_u32, 10_u32], 96_u32, 0.186_632, "p0"),
            ([64, 70], 134, 0.498_998, "p1"),
            ([164, 2], 166, 0.023_806_1, "p2"),
            ([314, 82], 396, 0.328_385, "the 200 individuals"),
        ] {
            assert_value(
                exp_het.of_var(&allele_counts(&of_each_allele), called_alleles, false),
                expected,
                OF_A_PRINTED_VALUE,
                &format!("the plain expected heterozygosity of var0000 in {pop}"),
            );
        }
    }

    /// No program outside the project prints the unbiased one, so what is
    /// checked is that the factor rebuilds it from plink2's own numbers: in
    /// p1 of `var0000` the called alleles are twice the 67 genotypes plink2
    /// counted, c = 134, and the plain 0.498998 times 134 / 133 is
    /// 0.502750, which is 4.5e-7 from pyNei's 0.5027494.
    #[test]
    fn the_unbiased_of_var0000_in_p1_is_the_plain_one_corrected() {
        let exp_het = ExpHet::new(2, 2, 20).unwrap();
        assert_value(
            exp_het.of_var(&allele_counts(&[64, 70]), 134, true),
            0.502_750,
            OF_A_PRINTED_VALUE,
            "the unbiased expected heterozygosity of var0000 in p1",
        );
    }

    /// The tetraploid variant of "What pyNei does that is odd": three
    /// individuals with the genotypes 0/1/2/3, 0/0/1/1 and 0/1/2/3, whose
    /// allele counts are 4, 4, 2 and 2 of 12 called alleles. pyNei corrects
    /// with the factor of a diploid population whatever the ploidy and
    /// gives 1.1685, above 1 and in no bin of the default histogram. popnei
    /// applies the correction of the exponent in hand and gives
    /// 1 - 48 / 11880, the products of four factors 24, 24, 0 and 0 over
    /// 12 · 11 · 10 · 9.
    #[test]
    fn the_tetraploid_variant_of_pynei_is_0_995960_unbiased() {
        let exp_het = ExpHet::new(4, 4, 1).unwrap();
        let counts = allele_counts(&[4, 4, 2, 2]);
        assert_value(
            exp_het.of_var(&counts, 12, true),
            0.995_960,
            OF_A_PRINTED_VALUE,
            "the unbiased expected heterozygosity of the tetraploid variant",
        );
    }

    /// The exponent is what the frequencies are raised to, and the ploidy
    /// is what turns the called alleles into called genotypes for the
    /// `min_num_individuals` test. The case is the one pyNei was measured
    /// on in "What pyNei does that is odd": four diploid individuals with
    /// the allele counts 5 and 3, asked for the exponent 4, whose plain
    /// value is 0.8276. Its 8 called alleles are 4 genotypes of the ploidy
    /// of the data and would be 2 of the exponent, so a threshold of 4
    /// keeps the value and one of 5 drops it.
    #[test]
    fn the_exponent_is_not_the_ploidy_that_counts_the_genotypes() {
        let counts = allele_counts(&[5, 3]);
        assert_value(
            ExpHet::new(4, 2, 1).unwrap().of_var(&counts, 8, false),
            0.8276,
            // The four digits the spec prints of pyNei's value.
            1e-4,
            "the plain expected heterozygosity at the exponent 4",
        );
        assert!(
            ExpHet::new(4, 2, 4)
                .unwrap()
                .of_var(&counts, 8, false)
                .is_some(),
            "8 called alleles are 4 genotypes of the ploidy 2"
        );
        assert_eq!(
            ExpHet::new(4, 2, 5).unwrap().of_var(&counts, 8, false),
            None
        );
    }

    /// At the exponent 1 the plain one is 1 minus the sum of the
    /// frequencies, which is 0 at every variant with a called allele, and
    /// so is the unbiased one, whose products have one factor each. pyNei
    /// computes the same.
    #[test]
    fn at_ploidy_1_a_variant_with_a_called_allele_is_0() {
        let exp_het = ExpHet::new(1, 1, 1).unwrap();
        let counts = allele_counts(&[3, 2]);
        for unbiased in [false, true] {
            assert_value(
                exp_het.of_var(&counts, 5, unbiased),
                0.0,
                OF_A_PRINTED_VALUE,
                &format!("the expected heterozygosity at ploidy 1, unbiased {unbiased}"),
            );
        }
    }

    /// A population that has called fewer than `min_num_individuals`
    /// genotypes has no value, and one with exactly that many keeps it. The
    /// called alleles are compared with that many times the ploidy, which
    /// is the same test pyNei makes by dividing.
    #[test]
    fn a_population_below_min_num_individuals_has_no_value() {
        let exp_het = ExpHet::new(2, 2, 20).unwrap();
        for unbiased in [false, true] {
            assert_eq!(
                exp_het.of_var(&allele_counts(&[29, 10]), 39, unbiased),
                None
            );
            assert!(
                exp_het
                    .of_var(&allele_counts(&[30, 10]), 40, unbiased)
                    .is_some(),
                "20 called genotypes, unbiased {unbiased}"
            );
        }
    }

    /// A population that has called nothing at a variant has no value,
    /// whatever the threshold. pyNei gives one for a whole chunk in which
    /// no individual is called anywhere, a plain 1 and an unbiased -0.0,
    /// because it takes the alleles it counts from the largest allele of
    /// the chunk; popnei works row by row and gives none.
    #[test]
    fn a_population_with_nothing_called_has_no_value() {
        let nothing_called = allele_counts(&[]);
        for unbiased in [false, true] {
            assert_eq!(
                ExpHet::new(2, 2, 0)
                    .unwrap()
                    .of_var(&nothing_called, 0, unbiased),
                None
            );
            assert_eq!(
                ExpHet::new(2, 2, 20)
                    .unwrap()
                    .of_var(&nothing_called, 0, unbiased),
                None
            );
        }
    }

    /// A diploid population with one called allele, which needs
    /// `min_num_individuals` at 0 to get this far, has fewer called alleles
    /// than the exponent: there are not two copies to draw. pyNei gives 0
    /// for the plain one and NaN for the unbiased one, and popnei the value
    /// and no value.
    #[test]
    fn one_called_allele_has_a_plain_value_and_no_unbiased_one() {
        let exp_het = ExpHet::new(2, 2, 0).unwrap();
        let counts = allele_counts(&[1]);
        assert_value(
            exp_het.of_var(&counts, 1, false),
            0.0,
            OF_A_PRINTED_VALUE,
            "the plain expected heterozygosity of one called allele",
        );
        assert_eq!(exp_het.of_var(&counts, 1, true), None);
    }

    /// An exponent of 0 gave, in the trial implementation the owner decided
    /// this with, a plain -2.0 and an unbiased NaN for the allele counts 2,
    /// 1 and 1, and one of 1e8 took 0.2 s for one variant; a ploidy of 0
    /// turns the `min_num_individuals` test off. The constructor refuses
    /// both, and each of the two arguments on its own.
    #[test]
    fn the_constructor_refuses_an_exponent_and_a_ploidy_of_0_and_of_256() {
        for out_of_range in [0, 256, usize::MAX] {
            let error = ExpHet::new(out_of_range, 2, 20).unwrap_err();
            assert!(
                matches!(&error, Error::StatPloidyOutOfRange { kind, value, largest }
                    if *kind == "exponent" && *value == out_of_range && *largest == 255),
                "the exponent {out_of_range}: {error:?}"
            );
            let error = ExpHet::new(2, out_of_range, 20).unwrap_err();
            assert!(
                matches!(&error, Error::StatPloidyOutOfRange { kind, value, largest }
                    if *kind == "ploidy" && *value == out_of_range && *largest == 255),
                "the ploidy {out_of_range}: {error:?}"
            );
        }
        assert!(ExpHet::new(1, 1, 20).is_ok());
        assert!(ExpHet::new(255, 255, 20).is_ok());
    }
}

#[cfg(test)]
mod distribs {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::{Path, PathBuf};

    use super::fixtures::{THE_SIX_VARIANTS, THE_THREE_VARIANTS};
    use super::{
        ExpHet, HistBins, Maf, ObsHet, PerVarDistribs, PerVarDistribsConfig, PerVarStat,
        PolyVarsStats, Pops, StatsDistrib, calc_per_var_distribs,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::filters::{FilteredReader, FilteringStats, VarFilter, VarFilteringCriterion};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, Needs};

    /// One unit of the last of the six digits the spec prints of a mean of
    /// the worked examples.
    const OF_A_PRINTED_MEAN: f64 = 1e-6;

    /// One unit of the last of the twelve digits the spec prints of the
    /// ratios of the polymorphism counts of the panel and of `many.vcf`.
    const OF_A_PRINTED_RATIO: f64 = 1e-12;

    /// A reader of the tests that gives the blocks it was built with, which
    /// is how the worked examples of the spec reach the pass: five diploid
    /// individuals named `i1` to `i5`, the individuals of both worked
    /// examples.
    #[derive(Debug)]
    struct GivenBlocks {
        individuals: Vec<String>,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
        /// What it was last asked to fill, which the test reads to see that
        /// the pass asks for the genotypes alone.
        needs: Needs,
    }

    impl GivenBlocks {
        /// The reader over `blocks`, which it gives in their order.
        fn of(blocks: Vec<Block>) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (1..=5).map(|number| format!("i{number}")).collect(),
                chroms,
                left,
                needs: Needs::ALL,
            }
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.left.pop())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            2
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }

    /// The blocks of `num_vars_per_block` variants that hold `variants`, the
    /// rows of a worked example, one after another.
    fn blocks_of(variants: &[[i8; 10]], num_vars_per_block: usize) -> Vec<Block> {
        variants
            .chunks(num_vars_per_block)
            .map(|of_the_block| Block {
                num_vars: of_the_block.len(),
                num_individuals: 5,
                ploidy: 2,
                gts: of_the_block.iter().flatten().copied().collect(),
                chrom: None,
                pos: None,
                id: None,
                alleles: None,
                qual: None,
            })
            .collect()
    }

    /// The reader over the six variants of the worked example of the pass,
    /// in blocks of `num_vars_per_block`.
    fn the_worked_example(num_vars_per_block: usize) -> GivenBlocks {
        GivenBlocks::of(blocks_of(&THE_SIX_VARIANTS, num_vars_per_block))
    }

    /// One population as `Pops::from_names` takes it.
    fn pop_of(name: &str, individuals: &[&str]) -> (String, Vec<String>) {
        (
            name.to_owned(),
            individuals
                .iter()
                .map(|individual| (*individual).to_owned())
                .collect(),
        )
    }

    /// The two populations of the worked examples, pop1 of i1 and i2 and
    /// pop2 of i3, i4 and i5, in that order.
    fn the_two_pops() -> Pops {
        let individuals: Vec<String> = (1..=5).map(|number| format!("i{number}")).collect();
        Pops::from_names(
            &[
                pop_of("pop1", &["i1", "i2"]),
                pop_of("pop2", &["i3", "i4", "i5"]),
            ],
            &individuals,
        )
        .expect("the two populations of the worked example")
    }

    /// The five statistics over `pops`, with `min_num_individuals`, the
    /// four bins of the worked examples and the default threshold of the
    /// polymorphism ratio.
    fn config_of(pops: Pops, min_num_individuals: u32) -> PerVarDistribsConfig {
        PerVarDistribsConfig {
            stats: vec![
                PerVarStat::ObsHet,
                PerVarStat::Maf,
                PerVarStat::ExpHet,
                PerVarStat::UnbiasedExpHet,
                PerVarStat::PolyVarsRatio,
            ],
            pops,
            bins: HistBins::linear(0.0, 1.0, 4).expect("the four bins of the worked example"),
            obs_het: ObsHet::new(min_num_individuals),
            maf: Maf::new(2, min_num_individuals).expect("the maf of diploid variants"),
            exp_het: ExpHet::new(2, 2, min_num_individuals)
                .expect("the expected heterozygosity of diploid variants"),
            poly_threshold: 0.95,
        }
    }

    /// That the mean of one population is the one the spec gives, within
    /// `tolerance`.
    fn assert_mean(distrib: Option<&StatsDistrib>, pop: usize, expected: f64, what: &str) {
        let distrib = distrib.unwrap_or_else(|| panic!("{what} was not calculated"));
        let Some(found) = distrib.mean(pop) else {
            panic!("{what} has no mean, and it is {expected}");
        };
        assert!(
            (found - expected).abs() <= OF_A_PRINTED_MEAN,
            "{what} is {found}, and it is {expected} within {OF_A_PRINTED_MEAN}"
        );
    }

    /// That the histogram of one population is the counts the spec gives,
    /// bin by bin.
    fn assert_hist(distrib: Option<&StatsDistrib>, pop: usize, expected: &[u64], what: &str) {
        let distrib = distrib.unwrap_or_else(|| panic!("{what} was not calculated"));
        assert_eq!(
            distrib.hist_counts(pop),
            expected,
            "the histogram of {what}"
        );
        let in_the_bins: u64 = expected.iter().copied().sum();
        assert!(
            in_the_bins <= distrib.num_vars_with_value(pop),
            "the histogram of {what} counts more variants than had a value"
        );
    }

    /// The means and the histograms of the worked example of the pass over
    /// its two populations, with `min_num_individuals` 1 and the four bins
    /// 0, 0.25, 0.5, 0.75 and 1, as the tables of "How it is verified" of
    /// the per variant distributions give them.
    ///
    /// The blocks of 6 and of 2 variants are the two the spec asks for: the
    /// pass adds the rows of every block into the same accumulators, so
    /// neither a count nor a mean depends on where the blocks were cut.
    #[test]
    fn the_means_and_the_histograms_of_the_worked_example_over_the_two_populations() {
        for num_vars_per_block in [6, 2] {
            let mut reader = the_worked_example(num_vars_per_block);
            let found = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
                .expect("the distributions of the worked example");

            assert_eq!(found.num_vars, 6);
            for (pop, means, hists) in [
                (
                    0,
                    [0.5, 0.6875, 0.375, 0.5],
                    [[1_u64, 0, 2, 1], [0, 1, 0, 3], [1, 2, 0, 1], [1, 0, 2, 1]],
                ),
                (
                    1,
                    [0.25, 0.729_167, 0.298_611, 0.383_333],
                    [[3, 0, 0, 1], [0, 1, 1, 2], [2, 1, 0, 1], [2, 0, 1, 1]],
                ),
            ] {
                for (index, (distrib, what)) in the_four_distribs(&found).iter().enumerate() {
                    let expected_mean = means.get(index).copied().expect("the mean of the table");
                    let expected_hist = hists.get(index).expect("the histogram of the table");
                    let what =
                        &format!("{what} of the population {pop}, blocks of {num_vars_per_block}");
                    assert_mean(*distrib, pop, expected_mean, what);
                    assert_hist(*distrib, pop, expected_hist, what);
                }
            }
        }
    }

    /// The four distributions of a result in the order of the tables of the
    /// spec: the observed heterozygosity, the maf, the plain expected
    /// heterozygosity and the unbiased one.
    fn the_four_distribs(found: &PerVarDistribs) -> [(Option<&StatsDistrib>, &'static str); 4] {
        [
            (found.obs_het.as_ref(), "the observed heterozygosity"),
            (found.maf.as_ref(), "the maf"),
            (found.exp_het.as_ref(), "the plain expected heterozygosity"),
            (
                found.unbiased_exp_het.as_ref(),
                "the unbiased expected heterozygosity",
            ),
        ]
    }

    /// The same example with no `pops`, which is the one population of the
    /// five individuals in the order of the source: the means and the
    /// histograms of the paragraph after the tables, over the variants 1,
    /// 2, 3 and 5.
    #[test]
    fn the_means_and_the_histograms_of_the_worked_example_with_no_pops() {
        for num_vars_per_block in [6, 2] {
            let mut reader = the_worked_example(num_vars_per_block);
            let found = calc_per_var_distribs(&mut reader, &config_of(Pops::all(5), 1))
                .expect("the distributions of the five individuals");

            let means = [0.395_833, 0.699_008, 0.378_107, 0.430_159];
            let hists = [[1_u64, 2, 0, 1], [0, 1, 0, 3], [2, 1, 0, 1], [1, 2, 0, 1]];
            for (index, (distrib, what)) in the_four_distribs(&found).iter().enumerate() {
                let what =
                    &format!("{what} of the five individuals, blocks of {num_vars_per_block}");
                assert_mean(
                    *distrib,
                    0,
                    means.get(index).copied().expect("the mean"),
                    what,
                );
                assert_hist(*distrib, 0, hists.get(index).expect("the histogram"), what);
                assert_eq!(
                    distrib.expect("the distribution").num_vars_with_value(0),
                    4,
                    "{what} is over the variants 1, 2, 3 and 5"
                );
            }
        }
    }

    /// The polymorphism counts and ratios of the worked example, from "How
    /// it is verified" of the polymorphism ratio: pop1 has the mafs 0.75,
    /// 0.75, 0.25 and 1, pop2 has 1, 1, 0.25 and 0.666667, and the five
    /// individuals together have 0.888889, 0.857143, 0.25 and 0.8.
    #[test]
    fn the_polymorphism_counts_and_ratios_of_the_worked_example() {
        for num_vars_per_block in [6, 2] {
            let mut reader = the_worked_example(num_vars_per_block);
            let found = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
                .expect("the distributions of the worked example");
            let poly = found
                .poly_vars_ratio
                .as_ref()
                .expect("the polymorphism counts");
            assert_counts(poly, 0, [3, 3, 4], "pop1");
            assert_ratios(poly, 0, 0.75, 1.0, "pop1");
            assert_counts(poly, 1, [2, 2, 4], "pop2");
            assert_ratios(poly, 1, 0.5, 1.0, "pop2");

            let mut reader = the_worked_example(num_vars_per_block);
            let found = calc_per_var_distribs(&mut reader, &config_of(Pops::all(5), 1))
                .expect("the distributions of the five individuals");
            let poly = found
                .poly_vars_ratio
                .as_ref()
                .expect("the polymorphism counts");
            assert_counts(poly, 0, [4, 4, 4], "the five individuals");
            assert_ratios(poly, 0, 1.0, 1.0, "the five individuals");
        }
    }

    /// The threshold is read at every variant: at 0.5 the only variant
    /// below it is the third, of maf 0.25 in both populations, so each has
    /// 1 polymorphic variant where the default 0.95 leaves 3 and 2, and
    /// neither the variable ones nor those with data change.
    #[test]
    fn a_poly_threshold_of_0_5_counts_the_variants_below_it_alone() {
        let mut reader = the_worked_example(6);
        let mut config = config_of(the_two_pops(), 1);
        config.poly_threshold = 0.5;
        let found = calc_per_var_distribs(&mut reader, &config)
            .expect("the distributions of the worked example");
        let poly = found
            .poly_vars_ratio
            .as_ref()
            .expect("the polymorphism counts");
        assert_counts(poly, 0, [1, 3, 4], "pop1");
        assert_ratios(poly, 0, 0.25, 0.333_333, "pop1");
        assert_counts(poly, 1, [1, 2, 4], "pop2");
        assert_ratios(poly, 1, 0.25, 0.5, "pop2");
    }

    /// That the three counts of one population are the ones the spec gives:
    /// the polymorphic variants, the variable ones and those with data.
    fn assert_counts(poly: &PolyVarsStats, pop: usize, expected: [u64; 3], what: &str) {
        assert_eq!(
            [
                poly.num_poly(pop),
                poly.num_variable(pop),
                poly.num_vars_with_data(pop)
            ],
            expected,
            "the polymorphic, the variable and the variants with data of {what}"
        );
    }

    /// That the two ratios of one population are the ones the spec gives,
    /// within one unit of the last digit it prints of them.
    fn assert_ratios(
        poly: &PolyVarsStats,
        pop: usize,
        over_the_data: f64,
        over_the_variables: f64,
        what: &str,
    ) {
        let found = poly.poly_ratio(pop).expect("the ratio over the data");
        assert!(
            (found - over_the_data).abs() <= OF_A_PRINTED_MEAN,
            "the ratio over the variants with data of {what} is {found}, and it is {over_the_data}"
        );
        let found = poly
            .poly_ratio_over_variables(pop)
            .expect("the ratio over the variable ones");
        assert!(
            (found - over_the_variables).abs() <= OF_A_PRINTED_MEAN,
            "the ratio over the variable variants of {what} is {found}, and it is \
             {over_the_variables}"
        );
    }

    /// The means of the worked example of the expected heterozygosity, its
    /// three variants over the same two populations with
    /// `min_num_individuals` 1: 0.3125 and 0.25 plain, 0.416667 and
    /// 0.333333 unbiased, each over the two variants that have a value.
    #[test]
    fn the_means_of_the_worked_example_of_the_expected_heterozygosity() {
        let mut reader = GivenBlocks::of(blocks_of(&THE_THREE_VARIANTS, 3));
        let found = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
            .expect("the distributions of the three variants");

        assert_eq!(found.num_vars, 3);
        assert_mean(found.exp_het.as_ref(), 0, 0.3125, "the plain one of pop1");
        assert_mean(found.exp_het.as_ref(), 1, 0.25, "the plain one of pop2");
        assert_mean(
            found.unbiased_exp_het.as_ref(),
            0,
            0.416_667,
            "the unbiased one of pop1",
        );
        assert_mean(
            found.unbiased_exp_het.as_ref(),
            1,
            0.333_333,
            "the unbiased one of pop2",
        );
    }

    /// The reference files of this module live at the root of the
    /// repository, beside the Python tests that read the same files.
    fn reference(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(name)
    }

    /// The reader over a VCF of the reference files, with every variant
    /// given and the blocks of the size popnei chose for its individuals.
    fn vcf_reader(name: &str) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: None,
        };
        VcfReader::<BufReader<File>>::from_path(&reference(name), options)
            .unwrap_or_else(|error| panic!("{name}: {error}"))
    }

    /// The populations `wanted` of a file of the reference, one line of an
    /// individual and its population, in the order `wanted` names them.
    ///
    /// The individuals of a population are not together in the file: the
    /// panel names its three in the order p0, p2, p1, and the population of
    /// a result is the one the caller asked for and not the one that
    /// happened to be named first.
    fn pops_of_the_file(name: &str, header: bool, individuals: &[String], wanted: &[&str]) -> Pops {
        let text = std::fs::read_to_string(reference(name)).expect("the populations of the file");
        let mut named: Vec<(String, Vec<String>)> = wanted
            .iter()
            .map(|pop| ((*pop).to_owned(), Vec::new()))
            .collect();
        for line in text.lines().skip(usize::from(header)) {
            let mut columns = line.split('\t');
            let (Some(individual), Some(pop)) = (columns.next(), columns.next()) else {
                panic!("the line `{line}` of {name} is not an individual and a population");
            };
            if let Some((_, of_the_pop)) = named.iter_mut().find(|(name, _)| name == pop) {
                of_the_pop.push(individual.to_owned());
            }
        }
        Pops::from_names(&named, individuals).expect("the populations of the file")
    }

    /// The nine polymorphism counts of the panel, from "How it is verified"
    /// of the polymorphism ratio: in `p0`, `p1` and `p2`, the variants
    /// whose maf is below 0.95, below 1, and that have one, with
    /// `min_num_individuals` 20. They come from the `--freq` reports of
    /// plink2 v2.0.0-a.7.7 and pyNei's `_calc_num_poly_vars` gives the
    /// same, and the two ratios of `p0` are those the spec prints to twelve
    /// digits.
    #[test]
    fn the_nine_polymorphism_counts_of_the_panel() {
        let mut reader = vcf_reader("stats/panel.vcf.gz");
        let pops = pops_of_the_file(
            "stats/panel_pops.txt",
            true,
            reader.individuals(),
            &["p0", "p1", "p2"],
        );
        let found = calc_per_var_distribs(&mut reader, &config_of(pops, 20))
            .expect("the distributions of the panel");

        assert_eq!(found.num_vars, 1200);
        let poly = found
            .poly_vars_ratio
            .as_ref()
            .expect("the polymorphism counts");
        assert_counts(poly, 0, [1112, 1173, 1200], "p0");
        assert_counts(poly, 1, [1101, 1177, 1200], "p1");
        assert_counts(poly, 2, [1093, 1184, 1200], "p2");
        assert_printed_ratio(poly.poly_ratio(0), 0.926_666_666_667, "p0 over the data");
        assert_printed_ratio(
            poly.poly_ratio_over_variables(0),
            0.947_996_589_940,
            "p0 over the variable ones",
        );
    }

    /// The six polymorphism counts of `many.vcf`, from the same part of the
    /// spec: `popA` of the first 20 individuals and `popB` of the other 30,
    /// with `min_num_individuals` 5, counted from bcftools 1.24, and the
    /// two ratios of `popA`.
    #[test]
    fn the_six_polymorphism_counts_of_many_vcf() {
        let mut reader = vcf_reader("vcf/many.vcf");
        let pops = pops_of_the_file(
            "stats/many_pops.txt",
            false,
            reader.individuals(),
            &["popA", "popB"],
        );
        let found = calc_per_var_distribs(&mut reader, &config_of(pops, 5))
            .expect("the distributions of many.vcf");

        assert_eq!(found.num_vars, 500);
        let poly = found
            .poly_vars_ratio
            .as_ref()
            .expect("the polymorphism counts");
        assert_counts(poly, 0, [477, 492, 500], "popA");
        assert_counts(poly, 1, [478, 493, 500], "popB");
        assert_printed_ratio(poly.poly_ratio(0), 0.954, "popA over the data");
        assert_printed_ratio(
            poly.poly_ratio_over_variables(0),
            0.969_512_195_122,
            "popA over the variable ones",
        );
    }

    /// That a ratio is the number the spec prints to twelve digits.
    fn assert_printed_ratio(found: Option<f64>, expected: f64, what: &str) {
        let found = found.unwrap_or_else(|| panic!("the ratio of {what} has no value"));
        assert!(
            (found - expected).abs() <= OF_A_PRINTED_RATIO,
            "the ratio of {what} is {found}, and it is {expected} within {OF_A_PRINTED_RATIO}"
        );
    }

    /// A population of 15 individuals at a `min_num_individuals` of 20 has
    /// no value of any statistic at any variant, whatever the data: every
    /// mean is none, every histogram counts nothing and the three
    /// polymorphism counts are 0, with both ratios none. It is the default
    /// threshold on the panel, which the pytest test of the spec asserts
    /// through Python.
    #[test]
    fn a_population_of_15_has_no_value_at_a_min_num_individuals_of_20() {
        let mut reader = vcf_reader("stats/panel.vcf.gz");
        let of_the_first_15: Vec<&str> = reader
            .individuals()
            .iter()
            .take(15)
            .map(String::as_str)
            .collect();
        let pops = Pops::from_names(&[pop_of("fifteen", &of_the_first_15)], reader.individuals())
            .expect("the population of 15 individuals");
        let found = calc_per_var_distribs(&mut reader, &config_of(pops, 20))
            .expect("the distributions of the panel");

        assert_eq!(found.num_vars, 1200);
        for (distrib, what) in the_four_distribs(&found) {
            let distrib = distrib.unwrap_or_else(|| panic!("{what} was not calculated"));
            assert_eq!(distrib.mean(0), None, "{what} has a mean");
            assert_eq!(distrib.num_vars_with_value(0), 0, "{what} has a variant");
            assert_eq!(distrib.hist_counts(0), [0, 0, 0, 0], "{what} counts a bin");
        }
        let poly = found
            .poly_vars_ratio
            .as_ref()
            .expect("the polymorphism counts");
        assert_counts(poly, 0, [0, 0, 0], "the population of 15");
        assert_eq!(poly.poly_ratio(0), None);
        assert_eq!(poly.poly_ratio_over_variables(0), None);
    }

    /// A pass over a source that holds no variant is refused, and the
    /// message says that the source holds none: a mean over no variant says
    /// nothing about a dataset, and the user has to know which of the two
    /// happened.
    #[test]
    fn a_pass_over_a_source_with_no_variant_is_refused() {
        let mut reader = GivenBlocks::of(Vec::new());
        let error = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
            .expect_err("a pass with no variant");
        assert!(
            matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source, filters }
                if *num_vars_of_the_source == 0 && filters.is_empty()),
            "{error:?}"
        );
        let message = error.to_string();
        assert!(message.contains("its source holds none"), "{message}");
    }

    /// A pass whose steps kept no variant is refused with the counts of
    /// each filter, which say where the variants went: the two variants of
    /// the block are the ones of the worked example with every genotype
    /// missing, and a missing data filter at 0 keeps neither.
    #[test]
    fn a_pass_whose_filter_kept_no_variant_is_refused_with_the_counts_of_the_filter() {
        let of_the_two = [
            *THE_SIX_VARIANTS.get(3).expect("the variant 4"),
            *THE_SIX_VARIANTS.get(5).expect("the variant 6"),
        ];
        let source = GivenBlocks::of(blocks_of(&of_the_two, 2));
        let filter = VarFilter::new(VarFilteringCriterion::MaxMissingRate(0.0))
            .expect("the missing data filter");
        let mut reader = FilteredReader::new(source, filter).expect("the chain of the pass");
        let error = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
            .expect_err("a pass with no variant");

        assert!(
            matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source, filters }
                if *num_vars_of_the_source == 2 && filters.len() == 1),
            "{error:?}"
        );
        let message = error.to_string();
        assert!(message.contains("its source gave 2"), "{message}");
        assert!(
            message.contains("`missing_data` filter was given 2 and kept 0"),
            "{message}"
        );
    }

    /// The threshold of the polymorphism ratio is a number from 0 to 1,
    /// both included, because a major allele frequency is one: the pass
    /// refuses anything else before it reads a variant.
    #[test]
    fn a_poly_threshold_that_is_not_a_number_from_0_to_1_is_refused() {
        for value in [-0.5, 1.5, f64::NAN] {
            let mut reader = the_worked_example(6);
            let mut config = config_of(the_two_pops(), 1);
            config.poly_threshold = value;
            let error = calc_per_var_distribs(&mut reader, &config)
                .expect_err("a threshold outside 0 to 1");
            assert!(
                matches!(error, Error::PolyThresholdOutOfRange { value: found }
                    if found.to_bits() == value.to_bits()),
                "{error:?}"
            );
        }
        // Both ends are taken.
        for value in [0.0, 1.0] {
            let mut reader = the_worked_example(6);
            let mut config = config_of(the_two_pops(), 1);
            config.poly_threshold = value;
            calc_per_var_distribs(&mut reader, &config).expect("a threshold at an end");
        }
    }

    /// A statistic that was not asked for is `None` in the result and
    /// changes none of the others: asking for fewer is a saving of work.
    #[test]
    fn a_statistic_that_was_not_asked_for_is_none_and_changes_no_other() {
        let mut reader = the_worked_example(6);
        let mut config = config_of(the_two_pops(), 1);
        config.stats = vec![PerVarStat::Maf];
        let found =
            calc_per_var_distribs(&mut reader, &config).expect("the maf of the worked example");

        assert!(found.obs_het.is_none());
        assert!(found.exp_het.is_none());
        assert!(found.unbiased_exp_het.is_none());
        assert!(found.poly_vars_ratio.is_none());
        assert_mean(found.maf.as_ref(), 0, 0.6875, "the maf of pop1");
        assert_hist(found.maf.as_ref(), 0, &[0, 1, 0, 3], "the maf of pop1");
    }

    /// The rows of a block are read on the threads of the pool the caller
    /// is in, and the chunks are of a fixed number of rows and are added in
    /// the order of the block, so any number of threads adds the same
    /// values in the same order: the counts are the same to the number and
    /// the means to the bit. "What could go wrong" of the plan asks for
    /// 1e-12 relative, and the bits are what the code gives.
    ///
    /// `many.vcf` in blocks of 150 variants spans four blocks of three
    /// chunks or fewer, so the rows of one block are read on several
    /// threads and the blocks follow one another.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_numbers_are_the_same_to_the_bit_in_pools_of_one_and_of_more_threads() {
        let of_the_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                let options = VcfOptions {
                    ploidy: 2,
                    only_passed: false,
                    num_vars_per_block: Some(150),
                };
                let mut reader =
                    VcfReader::<BufReader<File>>::from_path(&reference("vcf/many.vcf"), options)
                        .expect("many.vcf");
                let pops = pops_of_the_file(
                    "stats/many_pops.txt",
                    false,
                    reader.individuals(),
                    &["popA", "popB"],
                );
                calc_per_var_distribs(&mut reader, &config_of(pops, 5))
                    .expect("the distributions of many.vcf")
            })
        };
        let on_one = of_the_pool(1);
        assert_eq!(on_one.num_vars, 500);
        for threads in [2, 3, 4, 8, 16] {
            let on_more = of_the_pool(threads);
            assert_eq!(on_one.num_vars, on_more.num_vars);
            for ((of_one, what), (of_more, _)) in the_four_distribs(&on_one)
                .into_iter()
                .zip(the_four_distribs(&on_more))
            {
                let of_one = of_one.expect("the distribution on one thread");
                let of_more = of_more.expect("the distribution on more threads");
                for pop in 0..of_one.num_pops() {
                    assert_eq!(
                        of_one.hist_counts(pop),
                        of_more.hist_counts(pop),
                        "the histogram of {what} of the population {pop}, on {threads} threads"
                    );
                    assert_eq!(
                        of_one.num_vars_with_value(pop),
                        of_more.num_vars_with_value(pop)
                    );
                    let (Some(mean_of_one), Some(mean_of_more)) =
                        (of_one.mean(pop), of_more.mean(pop))
                    else {
                        panic!("{what} of the population {pop} has no mean");
                    };
                    // The bits and not a tolerance: the chunks are of a
                    // fixed number of rows and are added in the order of the
                    // block, so the threads add the same values in the same
                    // order, and a mean that differed in its last bit would
                    // say that the order had changed.
                    assert_eq!(
                        mean_of_one.to_bits(),
                        mean_of_more.to_bits(),
                        "the mean of {what} of the population {pop} is {mean_of_one} on one \
                         thread and {mean_of_more} on {threads}"
                    );
                }
            }
            let of_one = on_one
                .poly_vars_ratio
                .as_ref()
                .expect("the counts on one thread");
            let of_more = on_more
                .poly_vars_ratio
                .as_ref()
                .expect("the counts on more threads");
            for pop in 0..of_one.num_pops() {
                assert_counts(
                    of_more,
                    pop,
                    [
                        of_one.num_poly(pop),
                        of_one.num_variable(pop),
                        of_one.num_vars_with_data(pop),
                    ],
                    &format!("the population on {threads} threads"),
                );
            }
        }
    }
}
