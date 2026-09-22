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

use crate::error::{Error, Result};
use crate::filters::resolve_individuals;
use crate::io::vcf::MAX_PLOIDY;
use crate::variant::{AlleleCounts, GtCounts};

/// The name of the one population of a calculation that was given no
/// populations, inherited from pyNei's `DEF_POP_NAME`.
pub const DEFAULT_POP_NAME: &str = "pop";

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
    /// When `num_bins` is 0, when `start` is not below `end`, and when
    /// either of the two is NaN or infinite, which leaves every edge
    /// between them NaN.
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
/// A `num_bins` of 0, a `start` that is not below `end`, and a `start` or
/// an `end` that is NaN or infinite.
fn check_the_range(start: f64, end: f64, num_bins: usize) -> Result<()> {
    if num_bins == 0 {
        return Err(Error::HistWithNoBin);
    }
    if !start.is_finite() || !end.is_finite() || start >= end {
        return Err(Error::HistRangeNotGoingUp { start, end });
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
    use super::HistBins;
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
