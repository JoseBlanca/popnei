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
    count_the_genotype,
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

    /// The bins of the kind a user named, [`LINEAR_BINS`] of equal width or
    /// [`LOGARITHMIC_BINS`] of equal ratio.
    ///
    /// A user writes the kind as a string in Python and in TypeScript, so
    /// the name is read here and not in each binding crate.
    ///
    /// # Errors
    ///
    /// A `kind` that is neither of the two names, with both of them, and
    /// what the constructor the name picks refuses.
    pub fn of_kind(kind: &str, start: f64, end: f64, num_bins: usize) -> Result<HistBins> {
        match kind {
            LINEAR_BINS => HistBins::linear(start, end, num_bins),
            LOGARITHMIC_BINS => HistBins::logarithmic(start, end, num_bins),
            of_no_kind => Err(Error::HistBinsOfAnUnknownKind {
                kind: of_no_kind.to_owned(),
            }),
        }
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

    /// The same with the exponent the caller asked for, or the ploidy of
    /// the variants when they asked for none.
    ///
    /// It is what `ploidy=None` means in the pass in Python and in
    /// TypeScript, where the argument is the exponent and takes the ploidy
    /// of the variants when it is not given: the default is here and not in
    /// each binding crate.
    ///
    /// # Errors
    ///
    /// Those of [`ExpHet::new`].
    pub fn of_the_exponent_asked_for(
        exponent: Option<usize>,
        ploidy: usize,
        min_num_individuals: u32,
    ) -> Result<ExpHet> {
        ExpHet::new(exponent.unwrap_or(ploidy), ploidy, min_num_individuals)
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
pub(crate) fn raised(value: f64, times: u32) -> f64 {
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
pub(crate) fn checked_ploidy(kind: &'static str, value: usize) -> Result<u32> {
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

impl PerVarStat {
    /// The name of each of the five statistics, in the order of the
    /// variants above, which is the order of the fields of a result.
    ///
    /// The names are what a Python and a TypeScript user writes in `stats`,
    /// and each one is the field of the result that holds that statistic.
    /// They are here and not in the binding crates so that a rename is one
    /// change and not four.
    pub const NAMES: [&'static str; 5] = [
        "obs_het",
        "maf",
        "exp_het",
        "unbiased_exp_het",
        "poly_vars_ratio",
    ];

    /// The name a user writes for this statistic, which is the field of the
    /// result that holds it.
    #[must_use]
    pub fn name(self) -> &'static str {
        let of_the_five = match self {
            PerVarStat::ObsHet => 0,
            PerVarStat::Maf => 1,
            PerVarStat::ExpHet => 2,
            PerVarStat::UnbiasedExpHet => 3,
            PerVarStat::PolyVarsRatio => 4,
        };
        // The five names are there, one for each variant of the enum.
        PerVarStat::NAMES.get(of_the_five).copied().unwrap_or("")
    }

    /// The statistic a user named.
    ///
    /// # Errors
    ///
    /// A name that is of none of the five, with the five names.
    pub fn of_name(name: &str) -> Result<PerVarStat> {
        // The names are in [`PerVarStat::NAMES`] alone, in the order of the
        // variants, so a name that is renamed is renamed in one place.
        match PerVarStat::NAMES.iter().position(|known| *known == name) {
            Some(0) => Ok(PerVarStat::ObsHet),
            Some(1) => Ok(PerVarStat::Maf),
            Some(2) => Ok(PerVarStat::ExpHet),
            Some(3) => Ok(PerVarStat::UnbiasedExpHet),
            Some(4) => Ok(PerVarStat::PolyVarsRatio),
            Some(_) | None => Err(Error::StatOfAnUnknownName {
                name: name.to_owned(),
            }),
        }
    }
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

/// How many rows of a block one chunk of a pass reads, here and in the pass
/// over the populations of [`pop_dists`](crate::pop_dists).
///
/// The rows of a block are added up chunk by chunk and the chunks are added
/// together in the order of the block, so the sum of a statistic does not
/// depend on how many threads read the block, which rayon's own `sum` would
/// make it: it joins the parts in an order it chooses at run time. The
/// number is fixed for the same reason, and 64 rows of 1000 diploid
/// individuals are 128000 genotypes, enough work for one task of rayon.
pub(crate) const ROWS_PER_CHUNK: usize = 64;

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
/// asked for them and gave none, a block of no variants and a block whose
/// individuals or ploidy are not the ones the reader says its source has,
/// each a defect of a reader too; and a pass that gave no variant, whether
/// its source holds none or its steps kept none of them.
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
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    while let Some(block) = reader.next_block()? {
        let alleles_per_var = alleles_per_var_of(&block, num_individuals, ploidy)?;
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

/// How many alleles one variant of a block holds, its individuals times its
/// ploidy, which is how the rows of the block are cut, after the checks the
/// passes of this module and the pass over the populations of
/// [`pop_dists`](crate::pop_dists) make of every block their reader gives
/// them.
///
/// `num_individuals` and `ploidy` are what the reader says its source has.
/// A pass reads the rows of every block as rows of one run over the
/// variants, so each block has to be of those individuals and of that
/// ploidy: a block of others is read one individual at the place of
/// another, or counted whole for a population of every individual of the
/// reader, and the numbers that come out say nothing about themselves.
///
/// # Errors
///
/// An array of the block that is not of the size the block states; a block
/// of no variant; a block that holds the genotypes of no individual,
/// because it has no individual or because its ploidy is 0, which is told
/// apart from the genotypes that nobody asked the reader for, since a block
/// is empty of them in the same way; a block the genotypes are not in,
/// which a pass asked its reader for; and a block of other individuals or
/// of another ploidy than the reader says its source has. Each of them is a
/// defect of the reader that gave the block.
pub(crate) fn alleles_per_var_of(
    block: &Block,
    num_individuals: usize,
    ploidy: usize,
) -> Result<usize> {
    // The rows are cut out of the genotypes by the sizes the block states,
    // so those sizes are checked before anything is read.
    block.check()?;
    if block.num_vars == 0 {
        return Err(Error::ReaderGaveABlockOfNoVariants);
    }
    let alleles_per_var = block.alleles_per_var()?;
    if alleles_per_var == 0 {
        return Err(Error::BlockWithNoGenotypeOfAVariant {
            num_individuals: block.num_individuals,
            ploidy: block.ploidy,
        });
    }
    if block.gts.is_empty() {
        return Err(Error::FieldsNotInTheBlock { fields: Needs::GTS });
    }
    if block.num_individuals != num_individuals || block.ploidy != ploidy {
        return Err(Error::BlocksDoNotFitTogether {
            num_individuals,
            ploidy,
            found_num_individuals: block.num_individuals,
            found_ploidy: block.ploidy,
        });
    }
    Ok(alleles_per_var)
}

/// How many alleles one chunk of a pass holds: [`ROWS_PER_CHUNK`] rows of
/// `alleles_per_var` alleles, and one allele at least, because a cut of 0
/// is what the standard library refuses with a panic.
pub(crate) fn alleles_of_a_chunk(alleles_per_var: usize) -> usize {
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
    // `count_alleles_of` clears before it counts, so the loop over the rows
    // allocates nothing for a variant. What the pass allocates is the
    // `Totals` of each chunk of rows: the `Vec` of its populations and, in
    // it, one `Vec` of bin counts for each population and each of the four
    // statistics that have bins, which for 50 populations is 200 of them
    // per chunk.
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

/// What a pass has counted of one individual: the variants at which its
/// genotype is missing and the ones at which it is heterozygous.
#[derive(Debug, Clone, Copy)]
struct OfAnIndividual {
    num_missing: u64,
    num_het: u64,
}

impl OfAnIndividual {
    /// The counts of an individual before any variant is read.
    fn none() -> OfAnIndividual {
        OfAnIndividual {
            num_missing: 0,
            num_het: 0,
        }
    }

    /// It counts one genotype of the individual, which
    /// [`count_the_genotype`] has read into `of_the_genotype`: a missing
    /// genotype, a heterozygous one, or neither.
    fn add_the_genotype(&mut self, of_the_genotype: GtCounts) {
        // One genotype of the pass, and a pass of more than
        // 18446744073709551615 variants reads more rows than any source
        // holds.
        self.num_missing = self
            .num_missing
            .saturating_add(u64::from(of_the_genotype.missing));
        self.num_het = self.num_het.saturating_add(u64::from(of_the_genotype.het));
    }

    /// It adds what one chunk of rows counted to what the pass has.
    fn add_the_chunk(&mut self, of_the_chunk: OfAnIndividual) {
        self.num_missing = self.num_missing.saturating_add(of_the_chunk.num_missing);
        self.num_het = self.num_het.saturating_add(of_the_chunk.num_het);
    }
}

/// What one pass of [`calc_per_individual_stats`] gives back: how many
/// variants it gave, and of each individual how many of them its genotype
/// is missing at and how many it is heterozygous at.
///
/// The individuals are in the order of the rows of the blocks, which is the
/// order of `individuals()` of the reader of the pass, and the binding
/// crates take their names from there.
#[derive(Debug, Clone)]
pub struct PerIndividualStats {
    /// The two counts of each individual, in the order of the rows.
    individuals: Vec<OfAnIndividual>,
    /// The variants the pass gave.
    num_vars: u64,
}

impl PerIndividualStats {
    /// How many individuals it holds the counts of, those of the pass.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.individuals.len()
    }

    /// The variants the pass gave, which the missing rate of an individual
    /// is over.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// The variants at which the genotype of one individual is missing, a
    /// half called genotype among them.
    ///
    /// `individual` is one of `0..num_individuals()`; a number at or beyond
    /// `num_individuals()` is no individual of this and has no count here.
    #[must_use]
    pub fn num_missing(&self, individual: usize) -> u64 {
        self.individuals
            .get(individual)
            .map_or(0, |counts| counts.num_missing)
    }

    /// The variants at which the genotype of one individual is called and
    /// its alleles are not all the same.
    ///
    /// `individual` is one of `0..num_individuals()`; a number at or beyond
    /// `num_individuals()` is no individual of this and has no count here.
    #[must_use]
    pub fn num_het(&self, individual: usize) -> u64 {
        self.individuals
            .get(individual)
            .map_or(0, |counts| counts.num_het)
    }

    /// The missing genotypes of one individual over the variants of the
    /// pass, the share of them at which it has no genotype.
    ///
    /// `individual` is one of `0..num_individuals()`; a number at or beyond
    /// `num_individuals()` is no individual of this and has no rate here,
    /// NaN. A rate of 0 there would read as an individual whose genotype
    /// was called at every variant, and this one is `f64` and not an
    /// `Option<f64>` as the heterozygosity rate is.
    #[must_use]
    pub fn missing_rate(&self, individual: usize) -> f64 {
        let Some(counts) = self.individuals.get(individual) else {
            return f64::NAN;
        };
        // A pass that gave no variant is an error, so a result of
        // `calc_per_individual_stats` holds one variant at least and this
        // never happens; what it keeps out is the NaN of 0 over 0.
        if self.num_vars == 0 {
            return 0.0;
        }
        // Every count of popnei is below 2^53, where a `f64` holds the
        // whole numbers exactly.
        counts.num_missing as f64 / self.num_vars as f64
    }

    /// The heterozygous genotypes of one individual over its called ones,
    /// the variants of the pass less the ones its genotype is missing at,
    /// and `None` when it has called none of them.
    ///
    /// pyNei's `calc_per_sample_stats` divides by every variant instead, so
    /// an individual with more missing data looks less heterozygous there;
    /// the owner decided on 22 September 2026 that popnei divides by the
    /// called genotypes, which is what plink2's `--het` gives.
    ///
    /// `individual` is one of `0..num_individuals()`; a number at or beyond
    /// `num_individuals()` is no individual of this and has no rate here.
    #[must_use]
    pub fn obs_het_rate(&self, individual: usize) -> Option<f64> {
        let counts = self.individuals.get(individual)?;
        // The missing genotypes of an individual are counted among the
        // variants of the pass, so the called ones are not below 0.
        let called = self.num_vars.checked_sub(counts.num_missing)?;
        if called == 0 {
            return None;
        }
        // Every count of popnei is below 2^53, where a `f64` holds the
        // whole numbers exactly.
        Some(counts.num_het as f64 / called as f64)
    }
}

/// The missing rate and the heterozygosity rate of every individual over
/// the variants `reader` gives, which is one pass over the source through
/// the steps the variants carry.
///
/// `reader` is the outermost reader of the chain of the pass, lent and not
/// taken, so that whoever built the chain reads the counts of its filters
/// from it when this returns: those counts and
/// [`PerIndividualStats::num_vars`] are the `pass_stats` of a result in
/// Python and in TypeScript. The pass asks the reader for the genotypes
/// alone.
///
/// The individuals are those of `individuals()` of the reader, in that
/// order, which is the order of the genotypes in the rows of its blocks.
///
/// # Errors
///
/// What the reader fails with; a block that holds no genotypes, which is a
/// reader that was asked for them and gave none, a block of no variants and
/// a block whose individuals or ploidy are not the ones the reader says its
/// source has, each a defect of a reader; a genotype with an allele below
/// the missing one; and a pass that gave no variant, whether its source
/// holds none or its steps kept none of them.
pub fn calc_per_individual_stats<R: BlockReader + ?Sized>(
    reader: &mut R,
) -> Result<PerIndividualStats> {
    // The two counts of an individual follow from its genotype at each
    // variant, so no column of a block is read and the reader is asked to
    // fill none of them.
    reader.set_needs(Needs::GTS);
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    // The two counts of every individual, which every chunk of rows of
    // every block is added into: what is kept from one block to the next
    // grows with the individuals and not with the variants.
    let mut counted = vec![OfAnIndividual::none(); num_individuals];
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        let alleles_per_var = alleles_per_var_of(&block, num_individuals, ploidy)?;
        count_the_block(&block, alleles_per_var, &mut counted)?;
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
    Ok(PerIndividualStats {
        individuals: counted,
        num_vars,
    })
}

/// It counts the genotypes of every row of a block into `counted`, one
/// entry for each individual of the block.
///
/// Natively the chunks of rows are read on the threads of rayon, as section
/// 3 of `docs/architecture.md` asks: no row reads another, each chunk
/// counts its own rows, and the chunks are added into `counted` in the
/// order of the block. The counts are integers and add up to the same
/// number in any order; the order is kept all the same, so that a chunk is
/// added once and once only wherever it was read. The threads are those of
/// the pool the caller is running in, and rayon's global pool only when the
/// caller is in none.
///
/// # Errors
///
/// A genotype with an allele below the missing one, which no reader of
/// popnei gives. The error is the one of the first row that has one,
/// wherever the threads found it: which of two bad rows a thread reaches
/// first depends on how the chunks were shared out, and a user who reports
/// a damaged file has to get the same message every time, so the rows are
/// read again, one after another, to find the first.
#[cfg(not(target_family = "wasm"))]
fn count_the_block(
    block: &Block,
    alleles_per_var: usize,
    counted: &mut [OfAnIndividual],
) -> Result<()> {
    use rayon::iter::ParallelIterator;
    use rayon::slice::ParallelSlice;

    let num_individuals = counted.len();
    let of_the_chunks: Result<Vec<Vec<OfAnIndividual>>> = block
        .gts
        .par_chunks(alleles_of_a_chunk(alleles_per_var))
        .map(|chunk| {
            // The one allocation of a chunk: two counts for each
            // individual, which grow neither with the rows of the chunk nor
            // with the variants of the pass.
            let mut of_the_chunk = vec![OfAnIndividual::none(); num_individuals];
            count_the_rows(chunk, alleles_per_var, block.ploidy, &mut of_the_chunk)?;
            Ok(of_the_chunk)
        })
        .collect();
    match of_the_chunks {
        Ok(of_the_chunks) => {
            for of_the_chunk in &of_the_chunks {
                add_the_chunk(counted, of_the_chunk);
            }
            Ok(())
        }
        // The second pass costs a read of the block, and it is made only
        // where the block is refused and nothing of it is given.
        Err(of_a_thread) => {
            let mut read_again = vec![OfAnIndividual::none(); num_individuals];
            match count_the_chunks_one_by_one(block, alleles_per_var, &mut read_again) {
                Err(of_the_first_row) => Err(of_the_first_row),
                // The rows are the same rows, so the second pass finds an
                // error too; the error of the threads is what is left if it
                // ever did not.
                Ok(()) => Err(of_a_thread),
            }
        }
    }
}

/// The same counts, with the chunks read one after another, which is what
/// wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn count_the_block(
    block: &Block,
    alleles_per_var: usize,
    counted: &mut [OfAnIndividual],
) -> Result<()> {
    count_the_chunks_one_by_one(block, alleles_per_var, counted)
}

/// The chunks of the block read one after another, each into counts of its
/// own that are added into `counted` before the next is read: what wasm
/// runs, and what the threads fall back on to find the first row that is an
/// error.
///
/// The chunks are the same chunks the threads read, and they are added in
/// the same order.
///
/// # Errors
///
/// Those of [`count_the_block`], at the first row that has one.
fn count_the_chunks_one_by_one(
    block: &Block,
    alleles_per_var: usize,
    counted: &mut [OfAnIndividual],
) -> Result<()> {
    let mut of_the_chunk = vec![OfAnIndividual::none(); counted.len()];
    for chunk in block.gts.chunks(alleles_of_a_chunk(alleles_per_var)) {
        of_the_chunk.fill(OfAnIndividual::none());
        count_the_rows(chunk, alleles_per_var, block.ploidy, &mut of_the_chunk)?;
        add_the_chunk(counted, &of_the_chunk);
    }
    Ok(())
}

/// It adds what one chunk of rows counted to what the pass has, one
/// individual at a time.
fn add_the_chunk(counted: &mut [OfAnIndividual], of_the_chunk: &[OfAnIndividual]) {
    for (of_the_pass, of_the_chunk) in counted.iter_mut().zip(of_the_chunk) {
        of_the_pass.add_the_chunk(*of_the_chunk);
    }
}

/// It counts the genotypes of every row of one chunk into `of_the_chunk`,
/// one entry for each individual.
///
/// `gts` holds whole rows of `alleles_per_var` alleles each, and
/// `alleles_per_var` is the individuals of the block times its ploidy,
/// which the pass checked is not 0, so neither cut of the rows is of 0
/// alleles. The genotypes of a row are in the order of the individuals, so
/// the genotype at a place in the row and the counts at that place in
/// `of_the_chunk` are of the same individual.
///
/// # Errors
///
/// A genotype with an allele below the missing one, at the first row of the
/// chunk that has one.
fn count_the_rows(
    gts: &[i8],
    alleles_per_var: usize,
    ploidy: usize,
    of_the_chunk: &mut [OfAnIndividual],
) -> Result<()> {
    for row in gts.chunks_exact(alleles_per_var) {
        for (genotype, of_the_individual) in row.chunks_exact(ploidy).zip(of_the_chunk.iter_mut()) {
            // What a missing and a heterozygous genotype are is written
            // once, in the `variant` module, and the counts of one variant
            // read it there too.
            let mut of_the_genotype = GtCounts::default();
            count_the_genotype(genotype, &mut of_the_genotype)?;
            of_the_individual.add_the_genotype(of_the_genotype);
        }
    }
    Ok(())
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

    /// A number at or beyond the populations is no population of the
    /// `Pops`: it has no name, no individual, and it is not every
    /// individual of the reader, which is what the three doc comments say.
    /// A caller walks `0..len()` and never asks for one, and what those
    /// three give beyond it is what keeps them out of a panic.
    #[test]
    fn a_number_beyond_the_populations_has_no_name_and_no_individual() {
        let of_every_individual = Pops::all(5);
        let named = Pops::from_names(
            &[pop_of("pop1", &["i1", "i2"]), pop_of("pop2", &["i3"])],
            &the_five_individuals(),
        )
        .unwrap();
        for (pops, beyond) in [(&of_every_individual, 1), (&named, 2), (&named, 9)] {
            assert_eq!(pops.name(beyond), "");
            assert_eq!(pops.individuals(beyond), [] as [usize; 0]);
            assert!(!pops.is_all(beyond));
        }
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
mod per_var_stat {
    use super::PerVarStat;
    use crate::error::Error;

    /// The name of each statistic is what a user writes in `stats` and the
    /// field of the result that holds it, and `of_name` gives back the
    /// statistic of each name. The names live in `NAMES` alone, so a rename
    /// is one change; the literals here are the names of the fields of
    /// `PerVarDistribs` in Python and in TypeScript.
    #[test]
    fn each_name_is_the_name_of_the_statistic_it_gives_back() {
        assert_eq!(
            PerVarStat::NAMES,
            [
                "obs_het",
                "maf",
                "exp_het",
                "unbiased_exp_het",
                "poly_vars_ratio"
            ]
        );
        for (name, statistic) in [
            ("obs_het", PerVarStat::ObsHet),
            ("maf", PerVarStat::Maf),
            ("exp_het", PerVarStat::ExpHet),
            ("unbiased_exp_het", PerVarStat::UnbiasedExpHet),
            ("poly_vars_ratio", PerVarStat::PolyVarsRatio),
        ] {
            assert_eq!(
                PerVarStat::of_name(name).unwrap_or_else(|error| panic!("{name}: {error}")),
                statistic
            );
            assert_eq!(statistic.name(), name);
        }
    }

    /// A name that is of none of the five is refused, with the five names:
    /// a user who writes one of them wrong has to read which they are.
    #[test]
    fn a_name_of_no_statistic_is_refused_with_the_five() {
        for name in ["obs_hets", "OBS_HET", "", "exp_het "] {
            let error = PerVarStat::of_name(name).unwrap_err();
            assert!(
                matches!(&error, Error::StatOfAnUnknownName { name: found } if found == name),
                "{name}: {error:?}"
            );
            let message = error.to_string();
            for of_the_five in PerVarStat::NAMES {
                assert!(message.contains(of_the_five), "{message}");
            }
        }
    }
}

#[cfg(test)]
mod fixtures {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::{Path, PathBuf};

    use crate::block::{Block, BlockReader};
    use crate::error::Result;
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{
        AlleleCounts, ChromTable, GtCounts, Needs, count_alleles_of, count_gts, count_gts_of,
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

    /// The reference files live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    pub(super) fn reference(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference")
            .join(name)
    }

    /// The reader over a VCF of the reference files, with every variant
    /// given and the blocks of the size popnei chose for its individuals.
    pub(super) fn vcf_reader(name: &str) -> VcfReader<BufReader<File>> {
        vcf_reader_of(name, None)
    }

    /// The same reader, with blocks of `num_vars_per_block` variants when
    /// the test asks for a size of its own.
    pub(super) fn vcf_reader_of(
        name: &str,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block,
        };
        VcfReader::<BufReader<File>>::from_path(&reference(name), options)
            .unwrap_or_else(|error| panic!("{name}: {error}"))
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
            VcfReader::<BufReader<File>>::from_path(&reference("vcf/many.vcf"), options)
                .expect("many.vcf");
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

    /// A reader of the tests that gives the blocks it was built with, which
    /// is how the worked examples of the spec reach the pass: five diploid
    /// individuals named `i1` to `i5`, the individuals of both worked
    /// examples, unless the test asks for others.
    #[derive(Debug)]
    pub(super) struct GivenBlocks {
        individuals: Vec<String>,
        /// The ploidy it says its source has, which is the ploidy of the
        /// blocks it gives.
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks it has not given yet, the next one last.
        left: Vec<Block>,
        /// What it was last asked to fill, which the test reads to see that
        /// the pass asks for the genotypes alone.
        needs: Needs,
    }

    impl GivenBlocks {
        /// The fields the pass last asked it to fill.
        pub(super) fn needs(&self) -> Needs {
            self.needs
        }

        /// The reader over `blocks`, which it gives in their order, of the
        /// five diploid individuals of the worked examples.
        pub(super) fn of(blocks: Vec<Block>) -> GivenBlocks {
            GivenBlocks::of_a_source_of(5, 2, blocks)
        }

        /// The same reader over a source of `num_individuals` individuals,
        /// named `i1` to `iN`, of the ploidy `ploidy`, which is what it says
        /// its source has: a test of a pass over a ploidy that is not 2, and
        /// a test of what a pass does with a block that is not of the
        /// individuals or of the ploidy of its reader.
        pub(super) fn of_a_source_of(
            num_individuals: usize,
            ploidy: usize,
            blocks: Vec<Block>,
        ) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: (1..=num_individuals)
                    .map(|number| format!("i{number}"))
                    .collect(),
                ploidy,
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
            self.ploidy
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
    pub(super) fn blocks_of(variants: &[[i8; 10]], num_vars_per_block: usize) -> Vec<Block> {
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
}

#[cfg(test)]
mod hist {
    use super::{HistBins, LINEAR_BINS, LOGARITHMIC_BINS, MAX_NUM_BINS};
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

    /// The last edge is the end of the range itself and not the start plus
    /// the bins times the width, which is what `numpy.linspace` gives: the
    /// width of three bins from 0 to 0.9 is 0.3, and 0 plus three times it
    /// is 0.8999999999999999, one unit in the last place below the 0.9
    /// numpy and popnei end at. A value of exactly 0.9 is in the last bin
    /// with the one and outside the range with the other, so the two count
    /// the variants at the end of the range differently.
    #[test]
    fn the_last_edge_is_the_end_of_the_range_and_not_the_widths_added_up() {
        let bins = HistBins::linear(0.0, 0.9, 3).unwrap();
        assert_eq!(bins.edges(), [0.0, 0.3, 0.6, 0.9]);
        assert_eq!(bins.bin_of(0.9), Some(2));
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

    /// The kind a user writes picks the constructor: `linear` builds the
    /// bins of equal width and `logarithmic` those of equal ratio, and the
    /// two give what the constructors give. A kind that is neither is
    /// refused with both names, and `lineal`, which is how pyNei spells the
    /// first one, is one of those: the owner decided on 22 September 2026
    /// that popnei spells it `linear` and refuses the Spanish word.
    #[test]
    fn the_kind_of_bins_a_user_names_picks_the_constructor() {
        assert_eq!(
            HistBins::of_kind(LINEAR_BINS, 0.0, 1.0, 4)
                .expect("the bins of equal width")
                .edges(),
            [0.0, 0.25, 0.5, 0.75, 1.0]
        );
        assert_eq!(
            HistBins::of_kind(LOGARITHMIC_BINS, 0.01, 100.0, 4)
                .expect("the bins of equal ratio")
                .num_bins(),
            4
        );
        for kind in ["lineal", "log", "", "LINEAR"] {
            let error = HistBins::of_kind(kind, 0.0, 1.0, 4).unwrap_err();
            assert!(
                matches!(&error, Error::HistBinsOfAnUnknownKind { kind: found }
                    if found == kind),
                "{kind}: {error:?}"
            );
            let message = error.to_string();
            assert!(message.contains(LINEAR_BINS), "{message}");
            assert!(message.contains(LOGARITHMIC_BINS), "{message}");
        }
        // The kind names a constructor and refuses nothing of its own: what
        // that constructor refuses travels out as it is.
        let error = HistBins::of_kind(LOGARITHMIC_BINS, 0.0, 1.0, 4).unwrap_err();
        assert!(
            matches!(&error, Error::HistLogRangeNotAboveZero { .. }),
            "{error:?}"
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

    /// The three variants of four individuals of `test_obs_het_stats` of
    /// pyNei's `test/test_gt_counts.py`, over the one population of the
    /// four: `0/0 1/1 ./. 0/1`, which is 1 heterozygous genotype of 3
    /// called, `./1 0/0 0/1 1/0`, 2 of 3, whose first genotype is half
    /// called and missing, and four missing genotypes, which have no value.
    /// The values are the 1/3, 2/3 and none that test computes.
    #[test]
    fn the_three_values_of_the_test_of_pynei() {
        let obs_het = ObsHet::new(1);
        let of_the_three: [[i8; 8]; 3] = [
            [0, 0, 1, 1, -1, -1, 0, 1],
            [-1, 1, 0, 0, 0, 1, 1, 0],
            [-1, -1, -1, -1, -1, -1, -1, -1],
        ];
        let expected = [Some(1.0 / 3.0), Some(2.0 / 3.0), None];
        for (var, gts) in of_the_three.iter().enumerate() {
            let found = obs_het.of_var(gt_counts(gts));
            match expected[var] {
                Some(value) => assert_value(
                    found,
                    value,
                    OF_A_PRINTED_VALUE,
                    &format!("the observed heterozygosity of the variant {var} of pyNei's test"),
                ),
                None => assert_eq!(found, None, "the variant {var} of pyNei's test has a value"),
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
    use super::fixtures::{
        POP1, POP2, THE_SIX_VARIANTS, THE_THREE_VARIANTS, allele_counts, allele_counts_of,
        assert_value,
    };
    use super::{DEFAULT_HIST_RANGE, DEFAULT_NUM_BINS, ExpHet, HistBins};
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

    /// A caller that asks for no exponent gets the ploidy of the variants,
    /// which is what `ploidy=None` of the pass means in Python and in
    /// TypeScript. On the tetraploid variant above, 4, 4, 2 and 2 of 12,
    /// the exponent 4 gives 0.995960 and the exponent 2 another number, so
    /// the case says which of the two was taken.
    #[test]
    fn the_exponent_a_caller_asks_for_none_of_is_the_ploidy_of_the_variants() {
        let counts = allele_counts(&[4, 4, 2, 2]);
        let of_no_exponent = ExpHet::of_the_exponent_asked_for(None, 4, 1).unwrap();
        assert_value(
            of_no_exponent.of_var(&counts, 12, true),
            0.995_960,
            OF_A_PRINTED_VALUE,
            "the tetraploid variant with no exponent asked for",
        );
        let of_the_exponent_2 = ExpHet::of_the_exponent_asked_for(Some(2), 4, 1).unwrap();
        assert_eq!(
            of_the_exponent_2.of_var(&counts, 12, true),
            ExpHet::new(2, 4, 1).unwrap().of_var(&counts, 12, true)
        );
        assert_ne!(
            of_the_exponent_2.of_var(&counts, 12, true),
            of_no_exponent.of_var(&counts, 12, true)
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

    /// The 0 of the exponent 1 is the 0 of exact arithmetic: the
    /// frequencies are rounded, so their sum is not always 1. A haploid
    /// variant of nine individuals with the alleles 0, 1, 1, 1, 1, 1, 2, 3
    /// and 4 gives -2.220446049250313e-16, which numpy computes too, and
    /// that value falls below the range of the default histogram, so the
    /// variant counts in the mean and in no bin. popnei does not round it
    /// to 0, which would count it in the first bin where pyNei counts it in
    /// none.
    ///
    /// The value is compared to the bit: the frequencies and their sum are
    /// two of the four operations, which are rounded the same on every
    /// machine popnei runs on.
    #[test]
    fn at_ploidy_1_a_variant_can_fall_a_little_below_0() {
        let exp_het = ExpHet::new(1, 1, 1).unwrap();
        let counts = allele_counts(&[1, 5, 1, 1, 1]);
        for unbiased in [false, true] {
            assert_eq!(
                exp_het.of_var(&counts, 9, unbiased),
                Some(-2.220_446_049_250_313e-16),
                "the haploid variant of nine individuals, unbiased {unbiased}"
            );
        }
        let bins = HistBins::linear(DEFAULT_HIST_RANGE.0, DEFAULT_HIST_RANGE.1, DEFAULT_NUM_BINS)
            .expect("the default histogram");
        assert_eq!(bins.bin_of(-2.220_446_049_250_313e-16), None);
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

    use super::fixtures::{
        GivenBlocks, THE_SIX_VARIANTS, THE_THREE_VARIANTS, blocks_of, reference, vcf_reader,
        vcf_reader_of,
    };
    use super::{
        ExpHet, HistBins, Maf, ObsHet, PerVarDistribs, PerVarDistribsConfig, PerVarStat,
        PolyVarsStats, Pops, StatsDistrib, calc_per_var_distribs,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::Error;
    use crate::filters::{FilteredReader, VarFilter, VarFilteringCriterion};
    use crate::io::vcf::VcfReader;
    use crate::variant::Needs;

    /// One unit of the last of the six digits the spec prints of a mean of
    /// the worked examples.
    const OF_A_PRINTED_MEAN: f64 = 1e-6;

    /// One unit of the last of the twelve digits the spec prints of the
    /// ratios of the polymorphism counts of the panel and of `many.vcf`.
    const OF_A_PRINTED_RATIO: f64 = 1e-12;

    /// What two sizes of block are allowed to differ by in a mean, relative
    /// to the mean, which "What pyNei asserts, and the size of the blocks"
    /// of the spec gives: where the blocks were cut decides which rows are
    /// added up together, and the addition of floats is not associative.
    const OF_TWO_BLOCK_SIZES: f64 = 1e-12;

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
        for num_vars_per_block in [6, 2, 1] {
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
        for num_vars_per_block in [6, 2, 1] {
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
        for num_vars_per_block in [6, 2, 1] {
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

    /// The populations of `many.vcf` that the tests over it use, `popA` of
    /// the first 20 individuals and `popB` of the other 30.
    fn the_pops_of_many_vcf(reader: &VcfReader<BufReader<File>>) -> Pops {
        pops_of_the_file(
            "stats/many_pops.txt",
            false,
            reader.individuals(),
            &["popA", "popB"],
        )
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

    /// A value outside the range of the bins counts in the mean and in no
    /// bin, which "In Python and in TypeScript" of the pass states. The
    /// observed heterozygosities of pop1 in the worked example are 0.5,
    /// 0.5, 1 and 0, so with two bins from 0 to 0.5 the 1 of variant 3 is
    /// above the range: the mean is the 0.5 of the table of the spec, over
    /// the four variants that have a value, and the bins count three of
    /// them. A pass that dropped such a value from the mean too would give
    /// 1/3 here.
    #[test]
    fn a_value_above_the_range_of_the_bins_counts_in_the_mean_and_in_no_bin() {
        let mut reader = the_worked_example(6);
        let mut config = config_of(the_two_pops(), 1);
        config.bins = HistBins::linear(0.0, 0.5, 2).expect("two bins from 0 to 0.5");
        let found = calc_per_var_distribs(&mut reader, &config)
            .expect("the distributions of the worked example");

        let what = "the observed heterozygosity of pop1 over two bins from 0 to 0.5";
        assert_mean(found.obs_het.as_ref(), 0, 0.5, what);
        assert_hist(found.obs_het.as_ref(), 0, &[1, 2], what);
        assert_eq!(
            found
                .obs_het
                .as_ref()
                .expect("the observed heterozygosity")
                .num_vars_with_value(0),
            4,
            "{what} is over four variants and counts three of them in a bin"
        );
    }

    /// A block that holds the genotypes of no individual is refused for
    /// what it is and not as genotypes nobody asked the reader for: its
    /// `gts` are empty either way, and the two have nothing to do with each
    /// other. A reader of popnei gives neither, and the VCF reader refuses
    /// a header with no individual, so both say that a reader has a defect
    /// and the message has to name the right one.
    #[test]
    fn a_block_that_holds_the_genotypes_of_no_individual_is_refused_for_that() {
        let of_no_individual = Block {
            num_vars: 2,
            num_individuals: 0,
            ploidy: 2,
            gts: Vec::new(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let mut reader = GivenBlocks::of(vec![of_no_individual]);
        let error = calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
            .expect_err("a block of no individual");
        assert!(
            matches!(
                &error,
                Error::BlockWithNoGenotypeOfAVariant {
                    num_individuals: 0,
                    ploidy: 2
                }
            ),
            "{error:?}"
        );
    }

    /// The pass reads a row of genotypes as it is for a population of every
    /// individual of the reader, and it checks the width of the row before
    /// it does: a `Pops` built against another reader would otherwise count
    /// individuals the population does not hold and give a mean of them. A
    /// `Pops::all(5)` over the blocks of a reader of three individuals is
    /// the error of the counts of one variant, which name the individual
    /// beyond the row.
    #[test]
    fn a_pops_of_more_individuals_than_the_block_holds_is_refused() {
        let of_three_individuals = Block {
            num_vars: 2,
            num_individuals: 3,
            ploidy: 2,
            gts: vec![0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1],
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let mut reader = GivenBlocks::of_a_source_of(3, 2, vec![of_three_individuals]);
        let error = calc_per_var_distribs(&mut reader, &config_of(Pops::all(5), 1))
            .expect_err("a population of five individuals over a block of three");
        assert!(
            matches!(
                &error,
                Error::IndividualBeyondTheVariant {
                    individual: 3,
                    num_individuals: 3
                }
            ),
            "{error:?}"
        );
    }

    /// A block of other individuals than the reader says its source has is
    /// refused, and so is one of another ploidy. The pass reads the rows of
    /// every block as rows of one run over the variants, and a population
    /// of every individual of the reader is read as the whole row: a block
    /// of seven individuals after one of five would put seven genotypes
    /// into the counts of a population of five, which is a frequency over
    /// more data than the population holds.
    #[test]
    fn a_block_of_other_individuals_than_the_reader_says_is_refused() {
        let of_seven_individuals = Block {
            num_vars: 1,
            num_individuals: 7,
            ploidy: 2,
            gts: vec![0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1],
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let mut blocks = blocks_of(&THE_SIX_VARIANTS, 6);
        blocks.push(of_seven_individuals);
        let mut reader = GivenBlocks::of(blocks);
        let error = calc_per_var_distribs(&mut reader, &config_of(Pops::all(5), 1))
            .expect_err("a block of seven individuals of a reader of five");

        assert!(
            matches!(
                &error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 5,
                    ploidy: 2,
                    found_num_individuals: 7,
                    found_ploidy: 2
                }
            ),
            "{error:?}"
        );
    }

    /// The pass asks its reader for the genotypes alone, which "How it
    /// runs" of the per variant distributions states: the five statistics
    /// follow from them, and a reader that filled the columns of the
    /// chromosome, the position, the id, the alleles and the quality would
    /// read and hold what nothing reads.
    #[test]
    fn the_pass_asks_its_reader_for_the_genotypes_alone() {
        let mut reader = the_worked_example(6);
        calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
            .expect("the distributions of the worked example");
        assert_eq!(reader.needs(), Needs::GTS);
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
                let mut reader = vcf_reader_of("vcf/many.vcf", Some(150));
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
            assert_the_same_numbers(
                &on_one,
                &of_the_pool(threads),
                Agree::ToTheBit,
                &format!("{threads} threads against one"),
            );
        }
    }

    /// The chunks of a block read one after another, which is what wasm
    /// runs, give what the same chunks give on the threads of rayon, to the
    /// bit: they are the same chunks of the same rows and they are added in
    /// the same order. Nothing else compares the two paths, so without this
    /// the one a browser runs is only compiled.
    ///
    /// The block is the 500 variants of `many.vcf`, which are seven chunks
    /// of 64 rows and one of 52.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_chunks_of_a_block_read_one_by_one_give_what_the_threads_give() {
        use super::{Asked, Totals, add_the_block, add_the_chunks_one_by_one, the_distribs};

        let mut reader = vcf_reader_of("vcf/many.vcf", Some(500));
        reader.set_needs(Needs::GTS);
        let config = config_of(the_pops_of_many_vcf(&reader), 5);
        let block = reader
            .next_block()
            .expect("the block of many.vcf")
            .expect("many.vcf has variants");
        assert_eq!(block.num_vars, 500);
        let alleles_per_var = block.alleles_per_var().expect("the alleles of one variant");
        let asked = Asked::of(&config.stats);
        let of_the_path = |on_the_threads: bool| {
            let mut totals = Totals::of(config.pops.len(), config.bins.num_bins(), asked);
            if on_the_threads {
                add_the_block(&block, alleles_per_var, &config, asked, &mut totals)
                    .expect("the rows of the block on the threads");
            } else {
                add_the_chunks_one_by_one(&block, alleles_per_var, &config, asked, &mut totals)
                    .expect("the chunks of the block one by one");
            }
            the_distribs(&totals, &config, asked, 500)
        };

        assert_the_same_numbers(
            &of_the_path(true),
            &of_the_path(false),
            Agree::ToTheBit,
            "the chunks read one by one against the chunks read on the threads",
        );
    }

    /// The size of the blocks changes no histogram count and no mean beyond
    /// 1e-12 relative, which "What pyNei asserts, and the size of the
    /// blocks" promises. The sizes are compared against each other, and not
    /// each against the six digits the spec prints, which a drift of 1e-8
    /// between two sizes would pass.
    ///
    /// The worked example is six variants in blocks of 6, of 2 and of 1,
    /// and `many.vcf` is 500 in blocks of 500, of 150 and of 7, which cut
    /// the chunks of 64 rows in different places: where the blocks are cut
    /// decides which rows are added up together, and the addition of floats
    /// is not associative.
    #[test]
    fn the_size_of_the_blocks_changes_no_count_and_no_mean_beyond_1e_12() {
        let of_the_size = |num_vars_per_block| {
            let mut reader = the_worked_example(num_vars_per_block);
            calc_per_var_distribs(&mut reader, &config_of(the_two_pops(), 1))
                .expect("the distributions of the worked example")
        };
        let of_six = of_the_size(6);
        for num_vars_per_block in [2, 1] {
            assert_the_same_numbers(
                &of_six,
                &of_the_size(num_vars_per_block),
                Agree::WhereTheBlocksWereCut,
                &format!("blocks of {num_vars_per_block} variants against blocks of 6"),
            );
        }

        let of_the_vcf_size = |num_vars_per_block| {
            let mut reader = vcf_reader_of("vcf/many.vcf", Some(num_vars_per_block));
            let pops = the_pops_of_many_vcf(&reader);
            calc_per_var_distribs(&mut reader, &config_of(pops, 5))
                .expect("the distributions of many.vcf")
        };
        let of_the_whole_file = of_the_vcf_size(500);
        assert_eq!(of_the_whole_file.num_vars, 500);
        for num_vars_per_block in [150, 7] {
            assert_the_same_numbers(
                &of_the_whole_file,
                &of_the_vcf_size(num_vars_per_block),
                Agree::WhereTheBlocksWereCut,
                &format!("many.vcf in blocks of {num_vars_per_block} against one block"),
            );
        }
    }

    /// How closely two runs of the same pass have to agree.
    #[derive(Debug, Clone, Copy)]
    enum Agree {
        /// To the bit, which two runs that add the same values in the same
        /// order give: the threads of a pool, and the chunks of one block
        /// read one by one or on those threads.
        #[cfg_attr(
            target_family = "wasm",
            expect(
                dead_code,
                reason = "the two tests that compare runs adding the same values in the \
                          same order are the ones of the thread pools and of the chunks \
                          read one by one against the threads, and rayon is a dependency \
                          of the targets that are not wasm"
            )
        )]
        ToTheBit,
        /// Within 1e-12 relative, which is what the spec asks of two sizes
        /// of block, since where the blocks were cut decides which rows are
        /// added up together.
        WhereTheBlocksWereCut,
    }

    /// That two runs of the same pass gave the same counts, to the number,
    /// and the same means, as closely as `agree` asks.
    fn assert_the_same_numbers(
        left: &PerVarDistribs,
        right: &PerVarDistribs,
        agree: Agree,
        what: &str,
    ) {
        assert_eq!(left.num_vars, right.num_vars, "the variants of {what}");
        for ((of_left, statistic), (of_right, _)) in the_four_distribs(left)
            .into_iter()
            .zip(the_four_distribs(right))
        {
            let of_left = of_left.unwrap_or_else(|| panic!("{statistic} of {what}"));
            let of_right = of_right.unwrap_or_else(|| panic!("{statistic} of {what}"));
            for pop in 0..of_left.num_pops() {
                let of_them = &format!("{statistic} of the population {pop}, {what}");
                assert_eq!(
                    of_left.hist_counts(pop),
                    of_right.hist_counts(pop),
                    "the histogram of {of_them}"
                );
                assert_eq!(
                    of_left.num_vars_with_value(pop),
                    of_right.num_vars_with_value(pop),
                    "the variants with a value of {of_them}"
                );
                let (Some(mean_of_left), Some(mean_of_right)) =
                    (of_left.mean(pop), of_right.mean(pop))
                else {
                    panic!("{of_them} has no mean");
                };
                match agree {
                    Agree::ToTheBit => assert_eq!(
                        mean_of_left.to_bits(),
                        mean_of_right.to_bits(),
                        "the mean of {of_them} is {mean_of_left} and {mean_of_right}"
                    ),
                    Agree::WhereTheBlocksWereCut => assert!(
                        (mean_of_left - mean_of_right).abs()
                            <= OF_TWO_BLOCK_SIZES * mean_of_left.abs(),
                        "the mean of {of_them} is {mean_of_left} and {mean_of_right}"
                    ),
                }
            }
        }
        let of_left = left
            .poly_vars_ratio
            .as_ref()
            .expect("the counts of the one");
        let of_right = right
            .poly_vars_ratio
            .as_ref()
            .expect("the counts of the other");
        for pop in 0..of_left.num_pops() {
            assert_counts(
                of_right,
                pop,
                [
                    of_left.num_poly(pop),
                    of_left.num_variable(pop),
                    of_left.num_vars_with_data(pop),
                ],
                &format!("the population {pop}, {what}"),
            );
        }
    }
}

#[cfg(test)]
mod per_individual {
    use super::fixtures::{GivenBlocks, THE_SIX_VARIANTS, blocks_of, vcf_reader};
    use super::{PerIndividualStats, calc_per_individual_stats};
    use crate::block::Block;
    use crate::error::Error;
    use crate::variant::Needs;

    /// One unit of the last of the six digits the spec prints of a missing
    /// rate and of a heterozygosity rate.
    const OF_A_PRINTED_RATE: f64 = 1e-6;

    /// That one individual has the two counts and the two rates the spec
    /// gives it: the counts to the number, and the rates within the digits
    /// the spec prints of them.
    fn assert_the_numbers_of(
        found: &PerIndividualStats,
        individual: usize,
        counts: (u64, u64),
        rates: (f64, f64),
        what: &str,
    ) {
        let (num_missing, num_het) = counts;
        let (missing_rate, obs_het_rate) = rates;
        assert_eq!(
            found.num_missing(individual),
            num_missing,
            "the missing genotypes of {what}"
        );
        assert_eq!(
            found.num_het(individual),
            num_het,
            "the heterozygous genotypes of {what}"
        );
        let found_missing_rate = found.missing_rate(individual);
        assert!(
            (found_missing_rate - missing_rate).abs() <= OF_A_PRINTED_RATE,
            "the missing rate of {what} is {found_missing_rate}, and it is {missing_rate}"
        );
        let Some(found_obs_het_rate) = found.obs_het_rate(individual) else {
            panic!("{what} has no heterozygosity rate, and it is {obs_het_rate}");
        };
        assert!(
            (found_obs_het_rate - obs_het_rate).abs() <= OF_A_PRINTED_RATE,
            "the heterozygosity rate of {what} is {found_obs_het_rate}, and it is {obs_het_rate}"
        );
    }

    /// The name of one individual of the worked example, i1 to i5.
    fn individual_of_the_worked_example(individual: usize) -> String {
        format!("i{number}", number = individual.saturating_add(1))
    }

    /// A block as a reader gives one: `num_vars` variants of
    /// `num_individuals` individuals of the ploidy `ploidy`, whose
    /// genotypes are `gts`, and no column, which the pass asks for none of.
    fn block_of(num_vars: usize, num_individuals: usize, ploidy: usize, gts: Vec<i8>) -> Block {
        Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// The counts and the rates of the five individuals of the worked
    /// example, which "How it is verified" of the per individual statistics
    /// gives: the missing rates 2/6, 2/6, 2/6, 3/6 and 5/6 and the
    /// heterozygosity rates 1/4, 3/4, 1/4, 1/3 and 0/1. i5 has two half
    /// called genotypes, at the variants 1 and 2, and a half called
    /// genotype is missing and not heterozygous.
    ///
    /// The blocks of 6 and of 2 variants are the two the spec asks for, and
    /// blocks of 1 are a block for every row: the counts of each block are
    /// added into the same two counts of every individual, and the two
    /// divisions are made once at the end, so no rate depends on where the
    /// blocks were cut.
    #[test]
    fn the_counts_and_the_rates_of_the_five_individuals_of_the_worked_example() {
        for num_vars_per_block in [6, 2, 1] {
            let mut reader = GivenBlocks::of(blocks_of(&THE_SIX_VARIANTS, num_vars_per_block));
            let found = calc_per_individual_stats(&mut reader)
                .expect("the statistics of the worked example");

            assert_eq!(found.num_individuals(), 5);
            assert_eq!(found.num_vars(), 6);
            for (individual, counts, rates) in [
                (0, (2, 1), (0.333_333, 0.25)),
                (1, (2, 3), (0.333_333, 0.75)),
                (2, (2, 1), (0.333_333, 0.25)),
                (3, (3, 1), (0.5, 0.333_333)),
                (4, (5, 0), (0.833_333, 0.0)),
            ] {
                let what = format!(
                    "{name}, blocks of {num_vars_per_block}",
                    name = individual_of_the_worked_example(individual)
                );
                assert_the_numbers_of(&found, individual, counts, rates, &what);
            }
        }
    }

    /// The counts and the rates of `s000` and `s001` of the panel, the 1200
    /// biallelic diploid variants of 200 individuals of "How it is
    /// verified", read with the VCF reader: `s000` has 34 missing genotypes
    /// of 1200, 0.0283333, and 426 heterozygous of 1166 called, 0.365352;
    /// `s001` 44, 0.0366667, and 397 of 1156, 0.343426. They are the
    /// `MISSING_CT` of plink2's `--missing` and the `HET_CT` of its
    /// `--sample-counts`, over the `OBS_CT` of its `--het`.
    #[test]
    fn the_counts_and_the_rates_of_s000_and_s001_of_the_panel() {
        let mut reader = vcf_reader("stats/panel.vcf.gz");
        let found = calc_per_individual_stats(&mut reader).expect("the statistics of the panel");

        assert_eq!(found.num_individuals(), 200);
        assert_eq!(found.num_vars(), 1200);
        assert_the_numbers_of(&found, 0, (34, 426), (0.028_333, 0.365_352), "s000");
        assert_the_numbers_of(&found, 1, (44, 397), (0.036_667, 0.343_426), "s001");
    }

    /// The counts and the rates of `ind00` and `ind01` of `many.vcf`, the
    /// 500 variants of 50 diploid individuals with 257 half called
    /// genotypes and one variant in ten of three alleles: `ind00` has 29
    /// missing genotypes of 500, 0.058, and 201 heterozygous of 471 called,
    /// 0.426752; `ind01` 25, 0.05, and 195 of 475, 0.410526. They are the
    /// numbers of the same plink2 commands with `--vcf-half-call m`, and
    /// the `nMissing` and `nHets` of the `PSC` lines of `bcftools stats`.
    ///
    /// A heterozygous genotype of a third allele is heterozygous like any
    /// other, and a half called one is missing, as it is for pyNei.
    #[test]
    fn the_counts_and_the_rates_of_ind00_and_ind01_of_many_vcf() {
        let mut reader = vcf_reader("vcf/many.vcf");
        let found = calc_per_individual_stats(&mut reader).expect("the statistics of many.vcf");

        assert_eq!(found.num_individuals(), 50);
        assert_eq!(found.num_vars(), 500);
        assert_the_numbers_of(&found, 0, (29, 201), (0.058, 0.426_752), "ind00");
        assert_the_numbers_of(&found, 1, (25, 195), (0.05, 0.410_526), "ind01");
    }

    /// An individual with no called genotype has a missing rate of 1 and no
    /// heterozygosity rate, which is the NaN of the Python and the
    /// TypeScript results: 0 heterozygous genotypes of 0 called ones is no
    /// number.
    ///
    /// The two variants are the 3rd and the 6th of the worked example,
    /// `0/1 2/3 0/1 2/3 ./.` and `0/. ./. ./. ./. ./.`, so i1 to i4 are
    /// heterozygous at the first and missing at the second, a missing rate
    /// of 1/2 and a heterozygosity rate of 1/1, and i5 is missing at both.
    #[test]
    fn an_individual_with_no_called_genotype_has_a_missing_rate_of_1_and_no_het_rate() {
        let of_the_two = [
            *THE_SIX_VARIANTS.get(2).expect("the variant 3"),
            *THE_SIX_VARIANTS.get(5).expect("the variant 6"),
        ];
        let mut reader = GivenBlocks::of(blocks_of(&of_the_two, 2));
        let found =
            calc_per_individual_stats(&mut reader).expect("the statistics of the two variants");

        assert_eq!(found.num_vars(), 2);
        for individual in 0..4 {
            let what = individual_of_the_worked_example(individual);
            assert_the_numbers_of(&found, individual, (1, 1), (0.5, 1.0), &what);
        }
        assert_eq!(found.num_missing(4), 2);
        assert_eq!(found.num_het(4), 0);
        let missing_rate = found.missing_rate(4);
        assert!(
            (missing_rate - 1.0).abs() <= OF_A_PRINTED_RATE,
            "the missing rate of i5 is {missing_rate}, and it is 1"
        );
        assert_eq!(found.obs_het_rate(4), None, "the heterozygosity rate of i5");
    }

    /// A number at or beyond the individuals is no individual of the
    /// result: it has no missing genotype and no heterozygous one, 0 and 0,
    /// no heterozygosity rate, and a missing rate of NaN, which is what the
    /// five doc comments say. A caller walks `0..num_individuals()` and
    /// never asks for one, and what the five give beyond it is what keeps
    /// them out of a panic; a missing rate of 0.0 there would read as an
    /// individual whose genotype was called at every variant.
    #[test]
    fn a_number_beyond_the_individuals_has_no_count_and_no_rate() {
        let mut reader = GivenBlocks::of(blocks_of(&THE_SIX_VARIANTS, 6));
        let found =
            calc_per_individual_stats(&mut reader).expect("the statistics of the worked example");

        assert_eq!(found.num_individuals(), 5);
        for beyond in [5, 9, usize::MAX] {
            assert_eq!(found.num_missing(beyond), 0, "the {beyond}th individual");
            assert_eq!(found.num_het(beyond), 0, "the {beyond}th individual");
            let missing_rate = found.missing_rate(beyond);
            assert!(
                missing_rate.is_nan(),
                "the missing rate of the {beyond}th individual is {missing_rate}"
            );
            assert_eq!(
                found.obs_het_rate(beyond),
                None,
                "the {beyond}th individual"
            );
        }
    }

    /// The genotypes of a row are cut by the ploidy the reader says its
    /// source has, and not by 2: the rows of a tetraploid block hold four
    /// alleles of each individual.
    ///
    /// The two variants of the three tetraploid individuals are
    /// `0/0/0/0 0/0/1/1 ./././.` and `1/1/1/1 0/1/0/1 0/0/0/.`, so i1 is
    /// called and homozygous at both, i2 called and heterozygous at both,
    /// and i3 missing at both, the second of them a half called genotype.
    /// Cut by 2 the same rows would read the first six alleles as the
    /// genotypes of the three individuals, which gives i2 no heterozygous
    /// genotype and i3 one, and no missing genotype at all.
    #[test]
    fn the_rows_of_a_tetraploid_block_are_cut_by_the_ploidy_of_its_reader() {
        let of_three_tetraploid_individuals = block_of(
            2,
            3,
            4,
            vec![
                0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, -1,
            ],
        );
        let mut reader = GivenBlocks::of_a_source_of(3, 4, vec![of_three_tetraploid_individuals]);
        let found =
            calc_per_individual_stats(&mut reader).expect("the statistics of a tetraploid block");

        assert_eq!(found.num_individuals(), 3);
        assert_eq!(found.num_vars(), 2);
        assert_the_numbers_of(&found, 0, (0, 0), (0.0, 0.0), "i1 of the tetraploid block");
        assert_the_numbers_of(&found, 1, (0, 2), (0.0, 1.0), "i2 of the tetraploid block");
        assert_eq!(found.num_missing(2), 2, "the missing genotypes of i3");
        assert_eq!(found.num_het(2), 0, "the heterozygous genotypes of i3");
        let missing_rate = found.missing_rate(2);
        assert!(
            (missing_rate - 1.0).abs() <= OF_A_PRINTED_RATE,
            "the missing rate of i3 is {missing_rate}, and it is 1"
        );
        assert_eq!(found.obs_het_rate(2), None, "the heterozygosity rate of i3");
    }

    /// A pass over a source that holds no variant is refused, with the
    /// error of the pass of the per variant distributions: a rate over no
    /// variant says nothing about a dataset, and the message says whether
    /// the source held no variant or the steps kept none.
    #[test]
    fn a_pass_over_a_source_with_no_variant_is_refused() {
        let mut reader = GivenBlocks::of(Vec::new());
        let error = calc_per_individual_stats(&mut reader).expect_err("a pass with no variant");

        assert!(
            matches!(&error, Error::PassGaveNoVariant { num_vars_of_the_source, filters }
                if *num_vars_of_the_source == 0 && filters.is_empty()),
            "{error:?}"
        );
        let message = error.to_string();
        assert!(message.contains("its source holds none"), "{message}");
        // The message is the one of both passes, and a user of this one
        // reads it for a statistic of an individual and not of a variant.
        assert!(
            message.contains("a statistic of a pass is calculated over the variants it gives"),
            "{message}"
        );
    }

    /// A block of no variant is refused: every block a reader of popnei
    /// gives holds one variant at least, and a reader with no more variants
    /// gives no block, so a block of none says that the reader has a
    /// defect. Without the refusal the pass would read no row of it and go
    /// on to the next block with nothing said.
    #[test]
    fn a_block_of_no_variant_is_refused() {
        let mut reader = GivenBlocks::of(vec![block_of(0, 5, 2, Vec::new())]);
        let error = calc_per_individual_stats(&mut reader).expect_err("a block of no variant");

        assert!(
            matches!(&error, Error::ReaderGaveABlockOfNoVariants),
            "{error:?}"
        );
    }

    /// A block of other individuals than the reader says its source has is
    /// refused. The two counts of an individual are of the genotype at its
    /// place in every row of the pass, so a block of three individuals
    /// after one of five would count the individual 0 of the second block
    /// into the counts of the individual 0 of the first and leave the
    /// individuals 3 and 4 of the pass with the variants of the first block
    /// alone: three rates over one number of variants and two over
    /// another, with nothing to show it.
    #[test]
    fn a_block_of_other_individuals_than_the_reader_says_is_refused() {
        let of_three_individuals = block_of(2, 3, 2, vec![0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1]);
        let mut reader =
            GivenBlocks::of(vec![block_of(1, 5, 2, vec![0; 10]), of_three_individuals]);
        let error =
            calc_per_individual_stats(&mut reader).expect_err("a block of other individuals");

        assert!(
            matches!(
                &error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 5,
                    ploidy: 2,
                    found_num_individuals: 3,
                    found_ploidy: 2
                }
            ),
            "{error:?}"
        );
    }

    /// A block of another ploidy than the reader says its source has is
    /// refused for the same reason: the rows of a tetraploid block of five
    /// individuals are twice as wide, and cut by the ploidy of the reader
    /// they would give the genotype of one individual to another.
    #[test]
    fn a_block_of_another_ploidy_than_the_reader_says_is_refused() {
        let mut reader = GivenBlocks::of(vec![block_of(1, 5, 4, vec![0; 20])]);
        let error = calc_per_individual_stats(&mut reader).expect_err("a block of another ploidy");

        assert!(
            matches!(
                &error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 5,
                    ploidy: 2,
                    found_num_individuals: 5,
                    found_ploidy: 4
                }
            ),
            "{error:?}"
        );
    }

    /// A block that holds the genotypes of no individual is refused for
    /// what it is and not as genotypes nobody asked the reader for: its
    /// `gts` are empty either way, and the two have nothing to do with each
    /// other. The reader here says its source has no individual too, so
    /// what is wrong with the block is that a variant of it has no
    /// genotype.
    #[test]
    fn a_block_that_holds_the_genotypes_of_no_individual_is_refused_for_that() {
        let mut reader = GivenBlocks::of_a_source_of(0, 2, vec![block_of(2, 0, 2, Vec::new())]);
        let error = calc_per_individual_stats(&mut reader).expect_err("a block of no individual");

        assert!(
            matches!(
                &error,
                Error::BlockWithNoGenotypeOfAVariant {
                    num_individuals: 0,
                    ploidy: 2
                }
            ),
            "{error:?}"
        );
    }

    /// A block of variants whose genotypes are not there is refused: the
    /// pass asks its reader for the genotypes and reads nothing else, so a
    /// reader that gives a block without them has a defect. Without the
    /// refusal the pass would count no genotype of that block and divide by
    /// its variants all the same, which is a missing rate too low for every
    /// individual.
    #[test]
    fn a_block_whose_genotypes_are_not_there_is_refused() {
        let mut reader = GivenBlocks::of(vec![block_of(2, 5, 2, Vec::new())]);
        let error = calc_per_individual_stats(&mut reader).expect_err("a block with no genotypes");

        assert!(
            matches!(&error, Error::FieldsNotInTheBlock { fields } if *fields == Needs::GTS),
            "{error:?}"
        );
    }
    /// The pass reads the genotypes of a block and no other field, so it
    /// asks its reader for the genotypes alone: a reader that filled the
    /// chromosome, the position, the id, the alleles and the quality would
    /// read and hold what nothing reads.
    #[test]
    fn the_pass_asks_its_reader_for_the_genotypes_alone() {
        let mut reader = GivenBlocks::of(blocks_of(&THE_SIX_VARIANTS, 6));
        calc_per_individual_stats(&mut reader).expect("the statistics of the worked example");

        assert_eq!(reader.needs(), Needs::GTS);
    }

    /// The rows of a block are read on the threads of the pool the caller
    /// is in, and the chunks are of a fixed number of rows and are added in
    /// the order of the block, so any number of threads counts the same
    /// genotypes into the same counts: the counts agree to the number and
    /// the rates to the bit.
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
    fn the_numbers_are_the_same_to_the_bit_in_pools_of_one_and_of_four_threads() {
        use super::fixtures::vcf_reader_of;

        let of_the_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                let mut reader = vcf_reader_of("vcf/many.vcf", Some(150));
                calc_per_individual_stats(&mut reader).expect("the statistics of many.vcf")
            })
        };
        let on_one = of_the_pool(1);
        assert_eq!(on_one.num_individuals(), 50);
        assert_eq!(on_one.num_vars(), 500);

        assert_the_same_numbers(&on_one, &of_the_pool(4), "four threads against one");
    }

    /// The chunks of a block read one after another, which is what wasm
    /// runs, count what the same chunks count on the threads of rayon, to
    /// the bit: they are the same chunks of the same rows and they are
    /// added in the same order. Nothing else compares the two paths, so
    /// without this the one a browser runs is only compiled.
    ///
    /// The block is the 500 variants of `many.vcf`, which are seven chunks
    /// of 64 rows and one of 52.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_chunks_of_a_block_read_one_by_one_give_what_the_threads_give() {
        use super::fixtures::vcf_reader_of;
        use super::{OfAnIndividual, count_the_block, count_the_chunks_one_by_one};
        use crate::block::BlockReader;

        let mut reader = vcf_reader_of("vcf/many.vcf", Some(500));
        reader.set_needs(Needs::GTS);
        let block = reader
            .next_block()
            .expect("the block of many.vcf")
            .expect("many.vcf has variants");
        assert_eq!(block.num_vars, 500);
        let alleles_per_var = block.alleles_per_var().expect("the alleles of one variant");
        let of_the_path = |on_the_threads: bool| {
            let mut counted = vec![OfAnIndividual::none(); 50];
            if on_the_threads {
                count_the_block(&block, alleles_per_var, &mut counted)
                    .expect("the rows of the block on the threads");
            } else {
                count_the_chunks_one_by_one(&block, alleles_per_var, &mut counted)
                    .expect("the chunks of the block one by one");
            }
            PerIndividualStats {
                individuals: counted,
                num_vars: 500,
            }
        };

        assert_the_same_numbers(
            &of_the_path(true),
            &of_the_path(false),
            "the chunks read one by one against the chunks read on the threads",
        );
    }

    /// That two results hold the same counts, to the number, and the same
    /// rates, to the bit, for every individual.
    ///
    /// The counts are integers and add up to the same number in any order;
    /// the two divisions are made once at the end, over those counts, so
    /// the rates are the same bits and not merely the same number within a
    /// tolerance.
    #[cfg(not(target_family = "wasm"))]
    fn assert_the_same_numbers(left: &PerIndividualStats, right: &PerIndividualStats, what: &str) {
        assert_eq!(left.num_individuals(), right.num_individuals(), "{what}");
        assert_eq!(left.num_vars(), right.num_vars(), "{what}");
        for individual in 0..left.num_individuals() {
            let of_them = format!("the individual {individual}, {what}");
            assert_eq!(
                left.num_missing(individual),
                right.num_missing(individual),
                "the missing genotypes of {of_them}"
            );
            assert_eq!(
                left.num_het(individual),
                right.num_het(individual),
                "the heterozygous genotypes of {of_them}"
            );
            assert_eq!(
                left.missing_rate(individual).to_bits(),
                right.missing_rate(individual).to_bits(),
                "the missing rate of {of_them}"
            );
            assert_eq!(
                left.obs_het_rate(individual).map(f64::to_bits),
                right.obs_het_rate(individual).map(f64::to_bits),
                "the heterozygosity rate of {of_them}"
            );
        }
    }
}
