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
/// genotypes of that ploidy, and for an allele below [`MISSING_ALLELE`],
/// which no reader of popnei gives.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "each count is raised by one at most once for each genotype of `gts`, and \
              the genotypes are at most as many as its alleles, which were checked above \
              to be a number a u32 holds"
)]
pub fn count_gts(gts: &[i8], ploidy: usize) -> Result<GtCounts> {
    // `checked_rem` gives `None` for a ploidy of 0, which is the other
    // thing refused here and what `chunks_exact` below would panic at.
    if gts.len().checked_rem(ploidy) != Some(0) || u32::try_from(gts.len()).is_err() {
        return Err(Error::GtsNotWholeGenotypes {
            num_alleles: gts.len(),
            ploidy,
        });
    }
    let mut counts = GtCounts::default();
    for genotype in gts.chunks_exact(ploidy) {
        let mut alleles = genotype.iter().copied();
        // A chunk of `chunks_exact` holds the ploidy, which is 1 at
        // least, so every genotype has a first allele.
        let Some(first) = alleles.next() else {
            continue;
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
    }
    Ok(counts)
}

/// One count for each allele a genotype can hold, from 0 to
/// [`MAX_ALLELE`], which [`count_alleles`] fills.
pub type AlleleCounts = [u32; 128];

/// It adds to `counts[a]` how often the allele a was called in the
/// genotypes of one variant, and gives how many alleles it added, the
/// called alleles.
///
/// An allele is counted wherever it was called, in a half called genotype
/// too, which is what `_count_each_allele` of pyNei counts over a chunk.
/// The caller clears `counts` between two variants and hands the same
/// array over again, so that a pass over a block allocates nothing.
///
/// # Errors
///
/// For an allele below [`MISSING_ALLELE`], which no reader of popnei
/// gives.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the called alleles and each entry of `counts` are raised by one at most \
              once for each allele of `gts`, which were checked above to be a number a \
              u32 holds, and the caller clears `counts` between two variants"
)]
pub fn count_alleles(gts: &[i8], counts: &mut AlleleCounts) -> Result<u32> {
    // The counts of one variant are u32, so a variant of more alleles
    // than a u32 holds is refused instead of counted into a number that
    // wrapped. The counts of the alleles read one allele at a time, which
    // is the ploidy the error names.
    if u32::try_from(gts.len()).is_err() {
        return Err(Error::GtsNotWholeGenotypes {
            num_alleles: gts.len(),
            ploidy: 1,
        });
    }
    let mut called_alleles = 0_u32;
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
        called_alleles += 1;
    }
    Ok(called_alleles)
}

#[cfg(test)]
mod tests {
    use super::{AlleleCounts, ChromTable, GtCounts, Needs, count_alleles, count_gts};
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

    /// The counts are the caller's array, and what the function adds to it
    /// is one variant: a caller that reads a second variant into the same
    /// array without clearing it gets the two of them together.
    #[test]
    fn count_alleles_adds_to_the_counts_it_is_given() {
        let mut counts: AlleleCounts = [0; 128];
        assert_eq!(count_alleles(&THE_SIX_VARIANTS[0], &mut counts).unwrap(), 9);
        assert_eq!(
            count_alleles(&THE_SIX_VARIANTS[4], &mut counts).unwrap(),
            10
        );
        assert_eq!(counts[0], 16);
        assert_eq!(counts[1], 3);
        assert_eq!(counts[2], 0);

        counts = [0; 128];
        assert_eq!(count_alleles(&[], &mut counts).unwrap(), 0);
        assert_eq!(counts[0], 0);
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
}
