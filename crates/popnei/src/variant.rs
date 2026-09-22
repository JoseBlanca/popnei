//! What every other module of popnei says about a variant, one site of the
//! genome with the genotype of every individual at it: which of its fields
//! a consumer wants, the table that turns the name of a chromosome into a
//! number, the missing allele, and the view of one variant of a block.
//!
//! The variants flow in blocks, which [`crate::block`] holds, and a
//! calculation that works variant by variant walks the [`VariantRef`] of
//! the block it was given: a view into its arrays that allocates nothing.
//!
//! A consumer says with a [`Needs`] which fields it wants, and the reader
//! may skip the rest; a block says with
//! [`Block::fields`](crate::block::Block::fields) which ones it holds, and
//! a consumer that depends on a field it did not get fails with the error
//! that names it. `docs/specs/variant.md` has the design and section 1 of
//! `docs/architecture.md` the reasons for it.

use std::collections::HashMap;
use std::fmt;
use std::ops::{BitOr, BitOrAssign};

use crate::block::AllelesColumn;
use crate::error::{Error, Result};

/// An allele that was not called, `.` in a VCF.
pub const MISSING_ALLELE: i8 = -1;

/// The largest allele number a genotype can hold. 0 is the reference
/// allele and 1 and above are the alternative ones, in the order of the
/// VCF, so a variant has at most 128 alleles.
pub const MAX_ALLELE: i8 = i8::MAX;

/// The name of each field, for the messages. In the order of the bits.
const NAMES_OF_THE_NEEDS: [(Needs, &str); 5] = [
    (Needs::GTS, "gts"),
    (Needs::CHROM_POS, "chrom and pos"),
    (Needs::ID, "id"),
    (Needs::ALLELES, "alleles"),
    (Needs::QUAL, "qual"),
];

/// Which fields of the variants a consumer wants, or which ones a block
/// holds: a set of the five fields, with union, [`Needs::contains`] and
/// [`Needs::difference`].
///
/// A reader is asked for a set with
/// [`set_needs`](crate::block::BlockReader::set_needs) and may skip every
/// field that is not in it. Most calculations want the genotypes alone.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Needs(u8);

impl Needs {
    /// The genotypes, `gts`.
    pub const GTS: Needs = Needs(1);
    /// The chromosome and the position, which go together because no
    /// reader gives one without the other.
    pub const CHROM_POS: Needs = Needs(2);
    /// The id of the variant.
    pub const ID: Needs = Needs(4);
    /// The reference allele and the alternative ones.
    pub const ALLELES: Needs = Needs(8);
    /// The quality of the variant.
    pub const QUAL: Needs = Needs(16);
    /// The five fields, built from the five constants, so that a field
    /// added to this set later cannot be left out of it.
    pub const ALL: Needs = Needs::GTS
        .union(Needs::CHROM_POS)
        .union(Needs::ID)
        .union(Needs::ALLELES)
        .union(Needs::QUAL);

    /// No field at all.
    #[must_use]
    pub const fn empty() -> Needs {
        Needs(0)
    }

    /// Whether every field of `fields` is in this set. An empty `fields` is
    /// in every set.
    #[must_use]
    pub const fn contains(self, fields: Needs) -> bool {
        self.0 & fields.0 == fields.0
    }

    /// Whether this set holds no field.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }

    /// The fields of either set.
    #[must_use]
    pub const fn union(self, other: Needs) -> Needs {
        Needs(self.0 | other.0)
    }

    /// The fields of this set that are not in `other`. A consumer that
    /// depends on the fields it asked for gets from it, with the fields of
    /// the block it was given, the ones that are not there.
    #[must_use]
    pub const fn difference(self, other: Needs) -> Needs {
        Needs(self.0 & !other.0)
    }
}

impl BitOr for Needs {
    type Output = Needs;

    fn bitor(self, other: Needs) -> Needs {
        self.union(other)
    }
}

impl BitOrAssign for Needs {
    fn bitor_assign(&mut self, other: Needs) {
        *self = self.union(other);
    }
}

impl fmt::Display for Needs {
    /// The name of each field of the set between backticks, `` `gts`,
    /// `qual` ``, and `nothing` when it is empty. The backticks are what
    /// tells the reader of a message where one name ends, since one of the
    /// five is `chrom and pos`.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut written = false;
        for (field, name) in NAMES_OF_THE_NEEDS {
            if self.contains(field) {
                if written {
                    formatter.write_str(", ")?;
                }
                write!(formatter, "`{name}`")?;
                written = true;
            }
        }
        if !written {
            formatter.write_str("nothing")?;
        }
        Ok(())
    }
}

impl fmt::Debug for Needs {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "Needs({self})")
    }
}

/// The names of the chromosomes of one reader, each with its number.
///
/// A variant holds its chromosome as a number into this table, one `u32`
/// per variant instead of a string, and the reader keeps the table. The
/// numbers are given in the order in which the names first appear among
/// the variants the reader gives, so two passes over the same source give
/// the same numbers.
#[derive(Debug)]
pub struct ChromTable {
    names: Vec<String>,
    numbers: HashMap<String, u32>,
}

impl ChromTable {
    /// A table with no name in it.
    #[must_use]
    pub fn new() -> ChromTable {
        ChromTable {
            names: Vec::new(),
            numbers: HashMap::new(),
        }
    }

    /// The number of `name`, which is added when it is not there yet.
    ///
    /// A table holds at most `u32::MAX` names, which no genome reaches.
    /// Beyond that the number is `u32::MAX`, the name is not kept, and
    /// [`ChromTable::name`] gives `None` for it.
    pub fn intern(&mut self, name: &str) -> u32 {
        if let Some(number) = self.numbers.get(name) {
            return *number;
        }
        let Ok(number) = u32::try_from(self.names.len()) else {
            return u32::MAX;
        };
        if number == u32::MAX {
            return u32::MAX;
        }
        self.names.push(name.to_string());
        self.numbers.insert(name.to_string(), number);
        number
    }

    /// The name of `number`, or `None` when the table has no such number.
    #[must_use]
    pub fn name(&self, number: u32) -> Option<&str> {
        let index = usize::try_from(number).ok()?;
        self.names.get(index).map(String::as_str)
    }

    /// How many names the table holds.
    #[must_use]
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Whether the table holds no name.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
}

impl Default for ChromTable {
    fn default() -> ChromTable {
        ChromTable::new()
    }
}

/// One variant of a block: its genotypes and its other fields as they lie
/// in the arrays of that block.
///
/// A calculation that works variant by variant walks
/// [`Block::variants`](crate::block::Block::variants), and the row helpers
/// take this view. It allocates nothing and copies nothing: every field is
/// a number or a slice of the block it came from.
///
/// Every method but [`VariantRef::gts`] gives `None` when the block has no
/// such column, which is a column that nobody asked for or that the source
/// could not give.
#[derive(Debug, Clone, Copy)]
pub struct VariantRef<'a> {
    gts: &'a [i8],
    chrom: Option<u32>,
    pos: Option<u64>,
    id: Option<&'a str>,
    qual: Option<f32>,
    /// The alleles of the whole block, and which variant of them this is:
    /// the texts of one variant are not a slice of a column.
    alleles: Option<(&'a AllelesColumn, usize)>,
}

impl<'a> VariantRef<'a> {
    /// The view of one variant, which only the `block` module builds: the
    /// fields of a variant are read out of the columns of its block, and a
    /// view that another crate could build would not be a view of one.
    pub(crate) fn new(
        gts: &'a [i8],
        chrom: Option<u32>,
        pos: Option<u64>,
        id: Option<&'a str>,
        qual: Option<f32>,
        alleles: Option<(&'a AllelesColumn, usize)>,
    ) -> VariantRef<'a> {
        VariantRef {
            gts,
            chrom,
            pos,
            id,
            qual,
            alleles,
        }
    }

    /// The genotypes of the variant, num_individuals x ploidy alleles,
    /// individual after individual: the alleles of the individual i are
    /// `gts[i * ploidy .. (i + 1) * ploidy]`. 0 is the reference allele, 1
    /// up to [`MAX_ALLELE`] the alternative ones, and [`MISSING_ALLELE`] an
    /// allele that was not called.
    ///
    /// Empty when the block was built without the genotypes.
    #[must_use]
    pub fn gts(&self) -> &'a [i8] {
        self.gts
    }

    /// The number of the chromosome of the variant, in the [`ChromTable`]
    /// of the reader the block came from.
    #[must_use]
    pub fn chrom(&self) -> Option<u32> {
        self.chrom
    }

    /// The position of the variant, 1 based as in a VCF.
    #[must_use]
    pub fn pos(&self) -> Option<u64> {
        self.pos
    }

    /// The id of the variant, empty when the variant has none.
    #[must_use]
    pub fn id(&self) -> Option<&'a str> {
        self.id
    }

    /// The quality of the variant, phred scaled as the QUAL of a VCF:
    /// minus ten times the base ten logarithm of the probability that
    /// there is no variant at that site, so 30 is one in a thousand.
    ///
    /// It is NaN for a variant whose source gives no quality, which is
    /// what the column of a block holds for one, so a caller asks
    /// `is_nan` before it compares the quality or puts it in a sum: NaN
    /// travels through arithmetic and comes out at the end with nothing
    /// to say where it came from.
    #[must_use]
    pub fn qual(&self) -> Option<f32> {
        self.qual
    }

    /// How many alleles the variant has, the reference one among them.
    #[must_use]
    pub fn num_alleles(&self) -> Option<usize> {
        self.alleles.map(|(column, var)| column.num_alleles(var))
    }

    /// The text of one allele of the variant, as the source gave it: `A`,
    /// `<DEL>`, `*`. The allele 0 is the reference one.
    ///
    /// `None` for an allele the variant does not have, as for a block that
    /// holds no alleles.
    #[must_use]
    pub fn allele(&self, allele: usize) -> Option<&'a str> {
        let (column, var) = self.alleles?;
        let text = column.allele(var, allele);
        match text.is_empty() {
            true => None,
            false => Some(text),
        }
    }
}

/// The counts of the genotypes of one variant: how many were called, how
/// many are missing and how many are heterozygous.
///
/// [`count_gts`] works them out over the genotypes of one variant. A
/// genotype is missing when one of its alleles at least is
/// [`MISSING_ALLELE`], so a half called genotype, `0/.` in a VCF, is
/// missing and not called, and it is heterozygous when it is called and
/// its alleles are not all the same. `called` and `missing` add up to the
/// individuals of the variant, and `het` is at most `called`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GtCounts {
    /// Genotypes with no missing allele.
    pub called: u32,
    /// Genotypes with one missing allele at least, the half called among
    /// them.
    pub missing: u32,
    /// Called genotypes whose alleles are not all the same.
    pub het: u32,
}

/// How many of the genotypes of one variant are called, missing and
/// heterozygous.
///
/// `gts` is the genotypes of one variant, the alleles of one individual
/// after those of the individual before it, `ploidy` alleles each: a row
/// of the genotypes of a block, which [`VariantRef::gts`] gives. What a
/// missing and a heterozygous genotype are is in [`GtCounts`], and it is
/// what `_calc_gt_is_missing` and `_calc_gt_is_het` of pyNei compute as
/// masks.
///
/// # Errors
///
/// For a ploidy of 0 and for genotypes that are not a whole number of
/// genotypes of that ploidy, for a variant of more alleles than a count of
/// them holds, and for an allele below [`MISSING_ALLELE`], which no reader
/// of popnei gives.
pub fn count_gts(gts: &[i8], ploidy: usize) -> Result<GtCounts> {
    num_individuals_of(gts, ploidy)?;
    // The counts are `u32`, so a variant of more alleles than one holds is
    // refused instead of counted into a number that wrapped.
    if u32::try_from(gts.len()).is_err() {
        return Err(Error::MoreAllelesThanACountHolds {
            num_alleles: gts.len(),
        });
    }
    let mut counts = GtCounts::default();
    for genotype in gts.chunks_exact(ploidy) {
        count_the_genotype(genotype, &mut counts)?;
    }
    Ok(counts)
}

/// How many individuals `gts` holds the genotypes of, which is its alleles
/// divided by the ploidy.
///
/// # Errors
///
/// A ploidy of 0, and alleles that are not a whole number of genotypes of
/// the ploidy.
fn num_individuals_of(gts: &[i8], ploidy: usize) -> Result<usize> {
    // `checked_rem` and `checked_div` give `None` for a ploidy of 0, which
    // is the other thing refused here and what `chunks_exact` would panic
    // at.
    if gts.len().checked_rem(ploidy) != Some(0) {
        return Err(Error::GtsNotWholeGenotypes {
            num_alleles: gts.len(),
            ploidy,
        });
    }
    gts.len()
        .checked_div(ploidy)
        .ok_or(Error::GtsNotWholeGenotypes {
            num_alleles: gts.len(),
            ploidy,
        })
}

/// It counts one genotype into `counts`: the genotype of one individual at
/// one variant, the ploidy alleles of it.
///
/// The `stats` module counts it for one individual over the variants of a
/// pass, as [`count_gts`] counts it for one variant over the individuals,
/// so what a missing and a heterozygous genotype are is written here alone.
///
/// # Errors
///
/// An allele below [`MISSING_ALLELE`], which no reader of popnei gives.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "each count is raised by one at most once for each genotype counted, and \
              each caller checks first that the genotypes it counts are a number a u32 \
              holds"
)]
pub(crate) fn count_the_genotype(genotype: &[i8], counts: &mut GtCounts) -> Result<()> {
    let mut alleles = genotype.iter().copied();
    // A genotype holds the ploidy, which is 1 at least, so it has a first
    // allele.
    let Some(first) = alleles.next() else {
        return Ok(());
    };
    if first < MISSING_ALLELE {
        return Err(Error::AlleleBelowTheMissingOne { allele: first });
    }
    let mut missing = first == MISSING_ALLELE;
    let mut all_the_same = true;
    for allele in alleles {
        if allele < MISSING_ALLELE {
            return Err(Error::AlleleBelowTheMissingOne { allele });
        }
        missing |= allele == MISSING_ALLELE;
        all_the_same &= allele == first;
    }
    if missing {
        counts.missing += 1;
    } else {
        counts.called += 1;
        if !all_the_same {
            counts.het += 1;
        }
    }
    Ok(())
}

/// One count for each allele a genotype can hold, from 0 to
/// [`MAX_ALLELE`], which [`count_alleles`] fills.
pub type AlleleCounts = [u32; 128];

/// It writes into `counts[a]` how often the allele a was called in the
/// genotypes of one variant, and gives how many alleles it counted, the
/// called alleles.
///
/// An allele is counted wherever it was called, in a half called genotype
/// too, which is what `_count_each_allele` of pyNei counts over a chunk.
/// `counts` is cleared here before the variant is counted, so what it
/// holds afterwards is the counts of that variant and of no other: the
/// caller hands the same array over for every variant, which is what keeps
/// a pass over a block from allocating, and clears nothing itself. A
/// caller that had to clear it and forgot would get two variants added
/// together, and an entry already at the largest number a `u32` holds
/// would wrap with nothing to show it.
///
/// # Errors
///
/// For a variant of more alleles than a count of them holds, and for an
/// allele below [`MISSING_ALLELE`], which no reader of popnei gives.
pub fn count_alleles(gts: &[i8], counts: &mut AlleleCounts) -> Result<u32> {
    // The counts of one variant are u32, so a variant of more alleles
    // than a u32 holds is refused instead of counted into a number that
    // wrapped.
    if u32::try_from(gts.len()).is_err() {
        return Err(Error::MoreAllelesThanACountHolds {
            num_alleles: gts.len(),
        });
    }
    counts.fill(0);
    let mut called_alleles = 0_u32;
    count_the_alleles(gts, counts, &mut called_alleles)?;
    Ok(called_alleles)
}

/// It counts the alleles of `gts` into `counts`, which it does not clear,
/// and raises `called_alleles` by the ones that were called.
///
/// # Errors
///
/// An allele below [`MISSING_ALLELE`], which no reader of popnei gives.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the called alleles and each entry of `counts` are raised by one at most \
              once for each allele counted, and each caller checks first that the \
              alleles it counts are a number a u32 holds"
)]
fn count_the_alleles(
    gts: &[i8],
    counts: &mut AlleleCounts,
    called_alleles: &mut u32,
) -> Result<()> {
    for &allele in gts {
        if allele == MISSING_ALLELE {
            continue;
        }
        // An allele of 0 or more is at most `MAX_ALLELE`, which is the
        // largest an i8 holds, and `counts` has an entry for each one up
        // to it, so what the `else` catches is an allele below the
        // missing one.
        let Some(count) = usize::try_from(allele)
            .ok()
            .and_then(|entry| counts.get_mut(entry))
        else {
            return Err(Error::AlleleBelowTheMissingOne { allele });
        };
        *count += 1;
        *called_alleles += 1;
    }
    Ok(())
}

/// How many of the genotypes of one variant that belong to the individuals
/// of one population are called, missing and heterozygous.
///
/// `gts` is the genotypes of one variant, `ploidy` alleles for each
/// individual of the reader, and `individuals` the index of each individual
/// of the population among them, which [`Pops::individuals`] gives. It
/// counts the genotypes of those individuals and of no other, by the rules
/// of [`count_gts`], which a population of every individual in the order of
/// the reader is counted by instead: that one reads the row as it is.
///
/// [`Pops::individuals`]: crate::stats::Pops::individuals
///
/// # Errors
///
/// Those of [`count_gts`] over the genotypes it counts, and an individual
/// at or beyond the ones the variant holds the genotypes of, which is a
/// defect of popnei: the indices of a population are resolved before any
/// variant is read.
pub fn count_gts_of(gts: &[i8], ploidy: usize, individuals: &[usize]) -> Result<GtCounts> {
    let num_individuals = num_individuals_of(gts, ploidy)?;
    refuse_more_alleles_than_a_count_holds(individuals.len(), ploidy)?;
    let mut counts = GtCounts::default();
    for &individual in individuals {
        let genotype = genotype_of(gts, ploidy, individual, num_individuals)?;
        count_the_genotype(genotype, &mut counts)?;
    }
    Ok(counts)
}

/// It writes into `counts[a]` how often the allele a was called in the
/// genotypes of one variant that belong to the individuals of one
/// population, and gives how many alleles it counted, the called alleles of
/// the population.
///
/// `gts` is the genotypes of one variant, `ploidy` alleles for each
/// individual of the reader, and `individuals` the index of each individual
/// of the population among them, which [`Pops::individuals`] gives. It
/// counts the alleles of those individuals and of no other, by the rules of
/// [`count_alleles`], `counts` cleared here among them, which a population
/// of every individual in the order of the reader is counted by instead:
/// that one reads the row as it is.
///
/// [`Pops::individuals`]: crate::stats::Pops::individuals
///
/// # Errors
///
/// Those of [`count_alleles`] over the genotypes it counts, and an
/// individual at or beyond the ones the variant holds the genotypes of,
/// which is a defect of popnei: the indices of a population are resolved
/// before any variant is read.
pub fn count_alleles_of(
    gts: &[i8],
    ploidy: usize,
    individuals: &[usize],
    counts: &mut AlleleCounts,
) -> Result<u32> {
    let num_individuals = num_individuals_of(gts, ploidy)?;
    refuse_more_alleles_than_a_count_holds(individuals.len(), ploidy)?;
    counts.fill(0);
    let mut called_alleles = 0_u32;
    for &individual in individuals {
        let genotype = genotype_of(gts, ploidy, individual, num_individuals)?;
        count_the_alleles(genotype, counts, &mut called_alleles)?;
    }
    Ok(called_alleles)
}

/// The genotype of one individual at one variant: the ploidy alleles of
/// `gts` that are theirs.
///
/// # Errors
///
/// An individual at or beyond the `num_individuals` the variant holds the
/// genotypes of.
#[expect(
    clippy::unnecessary_lazy_evaluations,
    reason = "the error is built only when the individual is beyond the variant: with \
              `ok_or` rustc builds it before it knows whether it is needed, and the \
              destructor of `Error`, which is out of line because other cases of it own \
              strings, is then called on every lookup of every individual of every \
              population"
)]
fn genotype_of(
    gts: &[i8],
    ploidy: usize,
    individual: usize,
    num_individuals: usize,
) -> Result<&[i8]> {
    individual
        .checked_mul(ploidy)
        .and_then(|first| Some(first..first.checked_add(ploidy)?))
        .and_then(|of_the_individual| gts.get(of_the_individual))
        .ok_or_else(|| Error::IndividualBeyondTheVariant {
            individual,
            num_individuals,
        })
}

/// It refuses a population whose genotypes hold more alleles than a count
/// of the alleles of one variant holds, and with them more genotypes than a
/// count of the genotypes holds, since a genotype is one allele at least.
///
/// A population of a pass holds each of its individuals once, so no
/// population of a source reaches this; a caller that counts an individual
/// more than once can.
///
/// # Errors
///
/// `num_individuals` times `ploidy` above what a `u32` holds.
fn refuse_more_alleles_than_a_count_holds(num_individuals: usize, ploidy: usize) -> Result<()> {
    // What `saturating_mul` saturates at is above what a count holds, which
    // is what is refused here, so it hides nothing: it is the number of the
    // message that is the largest a `usize` holds and not the alleles.
    let num_alleles = num_individuals.saturating_mul(ploidy);
    if u32::try_from(num_alleles).is_err() {
        return Err(Error::MoreAllelesThanACountHolds { num_alleles });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        AlleleCounts, ChromTable, GtCounts, Needs, count_alleles, count_alleles_of, count_gts,
        count_gts_of,
    };
    use crate::error::Error;

    /// The six variants of five diploid individuals of the worked example
    /// of `docs/specs/filters.md`, which the table of "How it is verified"
    /// of the counts of one variant gives the counts of. `-1` is an allele
    /// that was not called, the `.` of a VCF.
    const THE_SIX_VARIANTS: [[i8; 10]; 6] = [
        // 0/0 0/1 0/0 0/0 0/.
        [0, 0, 0, 1, 0, 0, 0, 0, 0, -1],
        // 0/0 0/1 0/0 ./. 0/.
        [0, 0, 0, 1, 0, 0, -1, -1, 0, -1],
        // 0/1 2/3 0/1 2/3 ./.
        [0, 1, 2, 3, 0, 1, 2, 3, -1, -1],
        // ./. ./. ./. ./. ./.
        [-1; 10],
        // 0/0 0/0 0/0 0/0 1/1
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        // 0/. ./. ./. ./. ./.
        [0, -1, -1, -1, -1, -1, -1, -1, -1, -1],
    ];

    /// The alleles that were counted in `gts`, each with its count, and
    /// the called alleles. The alleles that were not called are left out,
    /// so that a test writes the counts as the table of the spec does.
    fn alleles_counted(gts: &[i8]) -> (Vec<(usize, u32)>, u32) {
        let mut counts: AlleleCounts = [0; 128];
        let called_alleles = count_alleles(gts, &mut counts).unwrap();
        let counted = counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(allele, count)| (allele, *count))
            .collect();
        (counted, called_alleles)
    }

    /// The counts of `_calc_gt_is_missing` and `_calc_gt_is_het` of pyNei
    /// at ef0ca6e on these genotypes, as the table of the spec has them. A
    /// half called genotype, the last one of the first variant, is missing
    /// and is not het.
    #[test]
    fn count_gts_of_the_six_variants_of_the_worked_example() {
        let counts = |variant: usize| count_gts(&THE_SIX_VARIANTS[variant], 2).unwrap();
        let expected = |called, missing, het| GtCounts {
            called,
            missing,
            het,
        };
        assert_eq!(counts(0), expected(4, 1, 1));
        assert_eq!(counts(1), expected(3, 2, 1));
        assert_eq!(counts(2), expected(4, 1, 4));
        assert_eq!(counts(3), expected(0, 5, 0));
        assert_eq!(counts(4), expected(5, 0, 0));
        assert_eq!(counts(5), expected(0, 5, 0));
    }

    /// The counts of `_count_alleles_per_var` of pyNei at ef0ca6e on these
    /// genotypes, as the table of the spec has them. The called allele of
    /// a half called genotype is counted, which is why the first variant
    /// has 9 called alleles of 10 and the last one has 1.
    #[test]
    fn count_alleles_of_the_six_variants_of_the_worked_example() {
        let counted = |variant: usize| alleles_counted(&THE_SIX_VARIANTS[variant]);
        assert_eq!(counted(0), (vec![(0, 8), (1, 1)], 9));
        assert_eq!(counted(1), (vec![(0, 6), (1, 1)], 7));
        assert_eq!(counted(2), (vec![(0, 2), (1, 2), (2, 2), (3, 2)], 8));
        assert_eq!(counted(3), (vec![], 0));
        assert_eq!(counted(4), (vec![(0, 8), (1, 2)], 10));
        assert_eq!(counted(5), (vec![(0, 1)], 1));
    }

    /// A genotype is heterozygous when its alleles are not all the same at
    /// any ploidy, and it is missing when one allele of it at least was
    /// not called, so 0/./0/0 is missing and not het. The counts are the
    /// ones the spec gives for these three genotypes.
    #[test]
    fn count_gts_of_the_three_tetraploid_genotypes() {
        let gts = [0, 0, 0, 1, 1, 1, 1, 1, 0, -1, 0, 0];
        assert_eq!(
            count_gts(&gts, 4).unwrap(),
            GtCounts {
                called: 2,
                missing: 1,
                het: 1,
            }
        );
        // The counts of the alleles are of the alleles and take no
        // ploidy: these three genotypes hold six 0 and five 1, the called
        // allele of the half called genotype among them. The spec's table
        // has the allele counts of the six diploid variants, and these
        // two numbers are counted off the genotypes above.
        assert_eq!(alleles_counted(&gts), (vec![(0, 6), (1, 5)], 11));
    }

    /// The genotypes of a variant are one genotype of the ploidy for each
    /// individual, so the counts refuse the alleles that are not that.
    /// They come from a reader with a defect, and the numbers of the
    /// message are what says which reader.
    #[test]
    fn count_gts_refuses_a_ploidy_of_0_and_genotypes_that_are_not_whole() {
        let error = count_gts(&[0, 0, 0, 1], 0).unwrap_err();
        assert!(
            matches!(
                error,
                Error::GtsNotWholeGenotypes {
                    num_alleles: 4,
                    ploidy: 0
                }
            ),
            "{error}"
        );

        let error = count_gts(&[0, 0, 0, 1, 0], 2).unwrap_err();
        assert!(
            matches!(
                error,
                Error::GtsNotWholeGenotypes {
                    num_alleles: 5,
                    ploidy: 2
                }
            ),
            "{error}"
        );

        // A variant of no individual is a whole number of genotypes, none,
        // and is counted and not refused.
        assert_eq!(count_gts(&[], 2).unwrap(), GtCounts::default());
        // The ploidy of the counts is the one they were given, and a
        // haploid genotype is called and never het.
        assert_eq!(
            count_gts(&[0, 1, -1], 1).unwrap(),
            GtCounts {
                called: 2,
                missing: 1,
                het: 0,
            }
        );
    }

    /// An allele below the missing one would be counted as a called
    /// allele, so the counts of the genotypes refuse it instead of giving
    /// a number that says nothing about it.
    #[test]
    fn count_gts_refuses_an_allele_below_the_missing_one() {
        let error = count_gts(&[0, 0, -2, 0], 2).unwrap_err();
        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
            "{error}"
        );

        let error = count_gts(&[i8::MIN, 0], 2).unwrap_err();
        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: i8::MIN }),
            "{error}"
        );
    }

    /// The counts of the alleles have one entry for each allele from 0 on,
    /// and an allele below the missing one has no entry to go into.
    #[test]
    fn count_alleles_refuses_an_allele_below_the_missing_one() {
        let mut counts: AlleleCounts = [0; 128];
        let error = count_alleles(&[0, 1, -2, -1], &mut counts).unwrap_err();
        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
            "{error}"
        );
    }

    /// The counts are the caller's array, handed over for one variant
    /// after another, and the function clears it before it counts: what it
    /// leaves is the counts of the variant it was given, whatever the
    /// array held, so a caller never adds two variants together and no
    /// entry of it can wrap.
    #[test]
    fn count_alleles_clears_the_counts_it_is_given() {
        let mut counts: AlleleCounts = [0; 128];
        assert_eq!(count_alleles(&THE_SIX_VARIANTS[0], &mut counts).unwrap(), 9);
        assert_eq!(counts[0], 8);
        assert_eq!(counts[1], 1);

        // The same array again, with the counts of the first variant in
        // it: the second variant has 8 zeros and 2 ones of its own.
        assert_eq!(
            count_alleles(&THE_SIX_VARIANTS[4], &mut counts).unwrap(),
            10
        );
        assert_eq!(counts[0], 8);
        assert_eq!(counts[1], 2);
        assert_eq!(counts[2], 0);

        // An entry at the largest number a count holds is cleared like any
        // other, where adding to it would wrap.
        counts[0] = u32::MAX;
        counts[3] = 7;
        assert_eq!(count_alleles(&THE_SIX_VARIANTS[2], &mut counts).unwrap(), 8);
        assert_eq!(counts[0], 2);
        assert_eq!(counts[3], 2);

        // A variant of no allele leaves the array empty and not what was
        // in it.
        assert_eq!(count_alleles(&[], &mut counts).unwrap(), 0);
        assert_eq!(counts[0], 0);
        assert_eq!(counts[3], 0);
    }

    #[test]
    fn a_set_of_needs_contains_the_fields_it_was_built_from_and_no_other() {
        let wanted = Needs::GTS | Needs::ID;
        assert!(wanted.contains(Needs::GTS));
        assert!(wanted.contains(Needs::ID));
        assert!(wanted.contains(Needs::GTS | Needs::ID));
        assert!(!wanted.contains(Needs::CHROM_POS));
        assert!(!wanted.contains(Needs::ALLELES));
        assert!(!wanted.contains(Needs::QUAL));
        assert!(!wanted.contains(Needs::ALL));

        assert!(Needs::ALL.contains(Needs::GTS | Needs::CHROM_POS | Needs::ID));
        assert!(Needs::ALL.contains(Needs::ALLELES | Needs::QUAL));
        assert!(Needs::ALL.contains(Needs::ALL));

        let nothing = Needs::empty();
        assert!(nothing.is_empty());
        assert!(!wanted.is_empty());
        assert!(!nothing.contains(Needs::GTS));
        assert!(nothing.contains(nothing));
        assert!(wanted.contains(nothing));

        let with_the_quality = wanted | Needs::QUAL;
        assert!(with_the_quality.contains(Needs::GTS | Needs::ID | Needs::QUAL));
        assert_eq!(with_the_quality.difference(wanted), Needs::QUAL);
        assert_eq!(Needs::ALL.difference(Needs::ALL), nothing);
        assert_eq!(wanted.difference(Needs::QUAL), wanted);

        let mut asked_for = Needs::GTS;
        asked_for |= Needs::ALLELES;
        assert_eq!(asked_for, Needs::GTS | Needs::ALLELES);
        assert_ne!(asked_for, Needs::GTS);
    }

    #[test]
    fn a_chrom_table_numbers_the_names_in_the_order_they_first_appear() {
        let mut chroms = ChromTable::new();
        assert!(chroms.is_empty());
        assert_eq!(chroms.len(), 0);
        assert_eq!(chroms.name(0), None);

        assert_eq!(chroms.intern("chr2"), 0);
        assert_eq!(chroms.intern("chr1"), 1);
        assert_eq!(chroms.intern("chr2"), 0);
        assert_eq!(chroms.intern("scaffold_7"), 2);
        assert_eq!(chroms.intern("chr1"), 1);

        assert_eq!(chroms.len(), 3);
        assert!(!chroms.is_empty());
        assert_eq!(chroms.name(0), Some("chr2"));
        assert_eq!(chroms.name(1), Some("chr1"));
        assert_eq!(chroms.name(2), Some("scaffold_7"));
        assert_eq!(chroms.name(3), None);
        assert_eq!(chroms.name(u32::MAX), None);
    }

    /// The individuals of the two populations of the worked example of
    /// `docs/specs/stats.md`, as the indices of `Pops::individuals`: pop1
    /// is i1 and i2, the first two individuals of the source, and pop2 is
    /// i3, i4 and i5, the other three.
    const POP1: [usize; 2] = [0, 1];
    const POP2: [usize; 3] = [2, 3, 4];

    /// The alleles that were counted in the genotypes of `individuals`,
    /// each with its count, and the called alleles of the population, as
    /// `alleles_counted` gives them over every individual.
    fn alleles_counted_of(
        gts: &[i8],
        ploidy: usize,
        individuals: &[usize],
    ) -> (Vec<(usize, u32)>, u32) {
        let mut counts: AlleleCounts = [0; 128];
        let called_alleles = count_alleles_of(gts, ploidy, individuals, &mut counts).unwrap();
        let counted = counts
            .iter()
            .enumerate()
            .filter(|(_, count)| **count > 0)
            .map(|(allele, count)| (allele, *count))
            .collect();
        (counted, called_alleles)
    }

    /// The counts of "How it is verified" of the counts of one variant over
    /// a population of `docs/specs/stats.md`, on variant 1 of the worked
    /// example, `0/0 0/1 0/0 0/0 0/.`. The half called genotype of i5 is
    /// missing and is not het, as it is over every individual.
    #[test]
    fn count_gts_of_variant_1_over_the_two_populations_of_the_worked_example() {
        assert_eq!(
            count_gts_of(&THE_SIX_VARIANTS[0], 2, &POP1).unwrap(),
            GtCounts {
                called: 2,
                missing: 0,
                het: 1,
            }
        );
        assert_eq!(
            count_gts_of(&THE_SIX_VARIANTS[0], 2, &POP2).unwrap(),
            GtCounts {
                called: 2,
                missing: 1,
                het: 0,
            }
        );
    }

    /// The counts of the same item of the spec: the called allele of the
    /// half called genotype of i5 is counted for pop2, which is why pop2
    /// has 5 called alleles of its 6.
    #[test]
    fn count_alleles_of_variant_1_over_the_two_populations_of_the_worked_example() {
        assert_eq!(
            alleles_counted_of(&THE_SIX_VARIANTS[0], 2, &POP1),
            (vec![(0, 3), (1, 1)], 4)
        );
        assert_eq!(
            alleles_counted_of(&THE_SIX_VARIANTS[0], 2, &POP2),
            (vec![(0, 5)], 5)
        );
    }

    /// A population of every individual in the order of the source is the
    /// whole row, so the two functions give what the counts over the row
    /// give. A caller that has such a population calls those instead,
    /// which read the row as it is.
    #[test]
    fn count_gts_of_every_individual_in_order_gives_what_count_gts_gives() {
        let every_individual = [0, 1, 2, 3, 4];
        for gts in &THE_SIX_VARIANTS {
            assert_eq!(
                count_gts_of(gts, 2, &every_individual).unwrap(),
                count_gts(gts, 2).unwrap()
            );
        }
    }

    /// As above, for the alleles: the same counts and the same called
    /// alleles.
    #[test]
    fn count_alleles_of_every_individual_in_order_gives_what_count_alleles_gives() {
        let every_individual = [0, 1, 2, 3, 4];
        for gts in &THE_SIX_VARIANTS {
            assert_eq!(
                alleles_counted_of(gts, 2, &every_individual),
                alleles_counted(gts)
            );
        }
    }

    /// An individual the variant has no genotype for is a defect of
    /// popnei, and the error names the index and how many individuals the
    /// variant holds, which is what says which of the two is wrong.
    #[test]
    fn count_gts_of_refuses_an_individual_beyond_the_variant() {
        let error = count_gts_of(&THE_SIX_VARIANTS[0], 2, &[0, 5]).unwrap_err();
        assert!(
            matches!(
                error,
                Error::IndividualBeyondTheVariant {
                    individual: 5,
                    num_individuals: 5
                }
            ),
            "{error:?}"
        );
    }

    /// As above, for the alleles.
    #[test]
    fn count_alleles_of_refuses_an_individual_beyond_the_variant() {
        let mut counts: AlleleCounts = [0; 128];
        let error = count_alleles_of(&THE_SIX_VARIANTS[0], 2, &[7], &mut counts).unwrap_err();
        assert!(
            matches!(
                error,
                Error::IndividualBeyondTheVariant {
                    individual: 7,
                    num_individuals: 5
                }
            ),
            "{error:?}"
        );
    }

    /// The rules of the genotypes are those of `count_gts`: a ploidy of 0,
    /// and genotypes that are not a whole number of genotypes of the
    /// ploidy, are refused before any individual is looked for.
    #[test]
    fn count_gts_of_refuses_a_ploidy_of_0_and_genotypes_that_are_not_whole() {
        let error = count_gts_of(&[0, 0, 0, 1], 0, &[0]).unwrap_err();
        assert!(
            matches!(
                error,
                Error::GtsNotWholeGenotypes {
                    num_alleles: 4,
                    ploidy: 0
                }
            ),
            "{error:?}"
        );

        let mut counts: AlleleCounts = [0; 128];
        let error = count_alleles_of(&[0, 0, 0, 1, 0], 2, &[0], &mut counts).unwrap_err();
        assert!(
            matches!(
                error,
                Error::GtsNotWholeGenotypes {
                    num_alleles: 5,
                    ploidy: 2
                }
            ),
            "{error:?}"
        );
    }

    /// An allele below the missing one is refused where it is counted, and
    /// the two functions count the individuals of the population alone: an
    /// allele of another individual is not read, as pyNei does not read it
    /// either, since it takes the genotypes of the population out of the
    /// chunk before it counts them.
    #[test]
    fn count_gts_of_and_count_alleles_of_refuse_an_allele_below_the_missing_one_of_the_pop() {
        let gts = [0, 0, -2, 1];
        let mut counts: AlleleCounts = [0; 128];
        let error = count_gts_of(&gts, 2, &[1]).unwrap_err();
        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
            "{error:?}"
        );
        let error = count_alleles_of(&gts, 2, &[1], &mut counts).unwrap_err();
        assert!(
            matches!(error, Error::AlleleBelowTheMissingOne { allele: -2 }),
            "{error:?}"
        );
        assert_eq!(
            count_gts_of(&gts, 2, &[0]).unwrap(),
            GtCounts {
                called: 1,
                missing: 0,
                het: 0,
            }
        );
        assert_eq!(alleles_counted_of(&gts, 2, &[0]), (vec![(0, 2)], 2));
    }

    /// The counts of the alleles are cleared before the variant is
    /// counted, as `count_alleles` clears them: what the array holds
    /// afterwards is the counts of that population at that variant, and
    /// the caller hands the same array over for every variant and every
    /// population.
    #[test]
    fn count_alleles_of_clears_the_counts_it_is_given() {
        let mut counts: AlleleCounts = [0; 128];
        assert_eq!(
            count_alleles_of(&THE_SIX_VARIANTS[2], 2, &POP2, &mut counts).unwrap(),
            4
        );
        assert_eq!(counts[0], 1);
        assert_eq!(counts[1], 1);
        assert_eq!(counts[2], 1);
        assert_eq!(counts[3], 1);

        assert_eq!(
            count_alleles_of(&THE_SIX_VARIANTS[4], 2, &POP1, &mut counts).unwrap(),
            4
        );
        assert_eq!(counts[0], 4);
        assert_eq!(counts[1], 0);
        assert_eq!(counts[2], 0);
        assert_eq!(counts[3], 0);
    }
}
