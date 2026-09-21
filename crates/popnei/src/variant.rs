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

#[cfg(test)]
mod tests {
    use super::{ChromTable, Needs};

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
