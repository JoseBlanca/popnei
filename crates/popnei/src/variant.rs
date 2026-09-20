//! The variant, one site of the genome with the genotype of every
//! individual at it, and the trait of anything that gives variants.
//!
//! A reader, the VCF reader, the vars file reader or a filter over another
//! reader, gives variants one at a time. The consumer, a calculation or a
//! writer, owns one [`Variant`], lends it to the reader again and again,
//! and the reader fills it with the next variant and says whether there was
//! one. The buffers inside the variant are allocated once and refilled, so
//! a million variants cost no allocation after the first few.
//!
//! A consumer says with a [`Needs`] which fields it wants, and the reader
//! may skip the rest; after each read the variant says in its `filled`
//! which fields it really holds. `docs/specs/variant.md` has the design and
//! section 1 of `docs/architecture.md` the reasons for it.

use std::collections::HashMap;
use std::fmt;
use std::ops::{BitOr, BitOrAssign};

use crate::error::Result;

/// An allele that was not called, `.` in a VCF.
pub const MISSING_ALLELE: i8 = -1;

/// The largest allele number a genotype can hold. 0 is the reference
/// allele and 1 and above are the alternative ones, in the order of the
/// VCF, so a variant has at most 128 alleles.
pub const MAX_ALLELE: i8 = i8::MAX;

/// The name of each field, for the messages. In the order of the bits.
const FIELD_NAMES: [(Needs, &str); 5] = [
    (Needs::GTS, "gts"),
    (Needs::CHROM_POS, "chrom and pos"),
    (Needs::ID, "id"),
    (Needs::ALLELES, "alleles"),
    (Needs::QUAL, "qual"),
];

/// Which fields of a [`Variant`] a consumer wants, or which ones a variant
/// holds: a set of the five fields, with union, [`Needs::contains`] and
/// [`Needs::difference`].
///
/// A reader is asked for a set with `set_needs` and may skip every field
/// that is not in it. Most calculations want the genotypes alone.
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

    /// No field at all, which is what `filled` of a variant that was just
    /// cleared holds.
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
    /// depends on the fields it asked for gets from it the ones the reader
    /// did not fill.
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
        for (field, name) in FIELD_NAMES {
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

/// One variant: one site of the genome with the genotype of every
/// individual at it.
///
/// The consumer owns it and lends it to a reader, which clears it and
/// fills it. Its fields are public because every reader writes them and
/// every consumer reads them. A field that is not in `filled` holds its
/// empty value and never what the previous variant left in it, so a
/// consumer that depends on a field checks `filled` first.
#[derive(Debug)]
pub struct Variant {
    /// A number of the [`ChromTable`] of the reader that filled this
    /// variant.
    pub chrom: u32,
    /// 1 based, as in a VCF.
    pub pos: u64,
    /// num_individuals x ploidy alleles, individual after individual: the
    /// alleles of individual i are `gts[i * ploidy .. (i + 1) * ploidy]`.
    /// 0 is the reference allele, 1 up to [`MAX_ALLELE`] the alternative
    /// ones in the order of the VCF, and [`MISSING_ALLELE`] an allele that
    /// was not called.
    pub gts: Vec<i8>,
    /// Empty when the source gives no id for the variant, `.` in a VCF,
    /// which is not the same as an id that was not asked for: that one has
    /// no `ID` in `filled`.
    pub id: String,
    /// The reference allele first, then the alternative ones.
    pub alleles: Vec<String>,
    /// The quality of the variant, the QUAL column of a VCF, phred
    /// scaled: minus ten times the base ten logarithm of the probability
    /// that there is no variant at this site, so 30 is one in a thousand.
    /// `None` when the source gives none.
    pub qual: Option<f32>,
    /// What the reader filled in the last read.
    pub filled: Needs,
}

impl Variant {
    /// A variant with every field at its empty value and nothing in
    /// `filled`.
    #[must_use]
    pub fn new() -> Variant {
        Variant {
            chrom: 0,
            pos: 0,
            gts: Vec::new(),
            id: String::new(),
            alleles: Vec::new(),
            qual: None,
            filled: Needs::empty(),
        }
    }

    /// Every field to its empty value and `filled` to nothing. It keeps
    /// the capacity of the buffers, which is what lets a reader refill a
    /// million variants without allocating.
    pub fn clear(&mut self) {
        self.clear_but_the_alleles();
        self.alleles.clear();
    }

    /// Every field but `alleles` to its empty value and `filled` to
    /// nothing. The alleles are left as they are, for the reader that
    /// writes over their strings instead of dropping them; it takes them
    /// out of the variant itself. [`Variant::clear`] is this followed by
    /// emptying `alleles`, so the list of the fields is written once and a
    /// field that is added later cannot be cleared by one and not the
    /// other.
    pub fn clear_but_the_alleles(&mut self) {
        self.chrom = 0;
        self.pos = 0;
        self.gts.clear();
        self.id.clear();
        self.qual = None;
        self.filled = Needs::empty();
    }
}

impl Default for Variant {
    fn default() -> Variant {
        Variant::new()
    }
}

/// Anything that gives variants one at a time: the VCF reader, the vars
/// file reader, a filter over another reader.
///
/// `read_variant` clears the variant it is lent, fills it and returns
/// true, or returns false when the source has no more variants; after a
/// false every later call returns false. An error ends the reader, and
/// what a call after an error returns is not defined.
///
/// The trait can be used as a boxed trait object, `Box<dyn
/// VariantReader>`, which is how the two binding crates hold their reader,
/// because neither a pyo3 class nor a wasm-bindgen class can be generic.
/// It asks for `Send` because a read ahead thread moves a reader into
/// another thread.
pub trait VariantReader: Send {
    /// Fills `var` with the next variant and returns true, or returns
    /// false when there are no more.
    ///
    /// # Errors
    ///
    /// When the source cannot be read or what it holds is malformed.
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool>;

    /// The names of the individuals, in the order of their genotypes in
    /// `gts`.
    fn individuals(&self) -> &[String];

    /// How many alleles the genotype of one individual holds.
    fn ploidy(&self) -> usize;

    /// The names of the chromosomes seen so far, each with its number.
    fn chroms(&self) -> &ChromTable;

    /// Which fields the reader is asked to fill. The rest may be skipped,
    /// and the change holds from the next read on. [`Needs::ALL`] until it
    /// is called.
    fn set_needs(&mut self, needs: Needs);
}

impl<R: VariantReader + ?Sized> VariantReader for Box<R> {
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
        (**self).read_variant(var)
    }

    fn individuals(&self) -> &[String] {
        (**self).individuals()
    }

    fn ploidy(&self) -> usize {
        (**self).ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        (**self).chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        (**self).set_needs(needs);
    }
}

#[cfg(test)]
mod tests {
    use super::{ChromTable, MAX_ALLELE, MISSING_ALLELE, Needs, Variant, VariantReader};
    use crate::error::Result;

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

    #[test]
    fn a_cleared_variant_is_empty_and_keeps_the_capacity_of_its_buffers() {
        let mut var = Variant::new();
        var.chrom = 4;
        var.pos = 1_000_003;
        var.gts
            .extend_from_slice(&[0, 1, MISSING_ALLELE, MAX_ALLELE]);
        var.id.push_str("rs4711");
        var.alleles.push("A".to_string());
        var.alleles.push("TTG".to_string());
        var.qual = Some(37.5);
        var.filled = Needs::ALL;

        let gts_capacity = var.gts.capacity();
        let id_capacity = var.id.capacity();
        let alleles_capacity = var.alleles.capacity();
        assert!(gts_capacity >= 4);
        assert!(id_capacity >= 6);
        assert!(alleles_capacity >= 2);

        var.clear();

        assert_eq!(var.chrom, 0);
        assert_eq!(var.pos, 0);
        assert!(var.gts.is_empty());
        assert!(var.id.is_empty());
        assert!(var.alleles.is_empty());
        assert_eq!(var.qual, None);
        assert_eq!(var.filled, Needs::empty());

        assert_eq!(var.gts.capacity(), gts_capacity);
        assert_eq!(var.id.capacity(), id_capacity);
        assert_eq!(var.alleles.capacity(), alleles_capacity);
    }

    /// The reader of a VCF writes the alleles over the strings the variant
    /// holds, so it needs the fields emptied and the alleles left alone.
    #[test]
    fn a_variant_cleared_but_the_alleles_keeps_them_and_empties_the_rest() {
        let mut var = Variant::new();
        var.chrom = 4;
        var.pos = 1_000_003;
        var.gts.extend_from_slice(&[0, MISSING_ALLELE]);
        var.id.push_str("rs4711");
        var.alleles.push("A".to_string());
        var.alleles.push("TTG".to_string());
        var.qual = Some(37.5);
        var.filled = Needs::ALL;

        var.clear_but_the_alleles();

        assert_eq!(var.chrom, 0);
        assert_eq!(var.pos, 0);
        assert!(var.gts.is_empty());
        assert!(var.id.is_empty());
        assert_eq!(var.qual, None);
        assert_eq!(var.filled, Needs::empty());
        assert_eq!(var.alleles, ["A", "TTG"]);
    }

    /// A reader of two variants of three individuals, written here to try
    /// the trait. It fills what it is asked for and nothing else, and it
    /// gives the numbers of its chromosomes in the order in which the
    /// names first appear.
    struct TwoVariants {
        individuals: Vec<String>,
        chroms: ChromTable,
        needs: Needs,
        /// The variants still to give, the last one first.
        left: Vec<(&'static str, u64, Vec<i8>, &'static str)>,
    }

    impl TwoVariants {
        fn new() -> TwoVariants {
            TwoVariants {
                individuals: vec![
                    "ind_1".to_string(),
                    "ind_2".to_string(),
                    "ind_3".to_string(),
                ],
                chroms: ChromTable::new(),
                needs: Needs::ALL,
                left: vec![
                    ("chr1", 24, vec![1, 1, 0, MISSING_ALLELE, 0, 0], "rs2"),
                    ("chr2", 11, vec![0, 0, 0, 1, MISSING_ALLELE, 1], "rs1"),
                ],
            }
        }
    }

    impl VariantReader for TwoVariants {
        fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
            var.clear();
            let Some((chrom, pos, gts, id)) = self.left.pop() else {
                return Ok(false);
            };
            var.chrom = self.chroms.intern(chrom);
            var.pos = pos;
            var.filled = Needs::CHROM_POS;
            if self.needs.contains(Needs::GTS) {
                var.gts.extend_from_slice(&gts);
                var.filled |= Needs::GTS;
            }
            if self.needs.contains(Needs::ID) {
                var.id.push_str(id);
                var.filled |= Needs::ID;
            }
            Ok(true)
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
    }

    #[test]
    fn two_variants_are_read_through_a_boxed_reader() {
        let mut reader: Box<dyn VariantReader> = Box::new(TwoVariants::new());
        let mut var = Variant::new();

        assert_eq!(reader.individuals(), ["ind_1", "ind_2", "ind_3"]);
        assert_eq!(reader.ploidy(), 2);
        assert_eq!(reader.chroms().len(), 0);

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.pos, 11);
        assert_eq!(var.gts, [0, 0, 0, 1, MISSING_ALLELE, 1]);
        assert_eq!(var.id, "rs1");
        assert_eq!(var.filled, Needs::GTS | Needs::CHROM_POS | Needs::ID);
        assert_eq!(reader.chroms().name(var.chrom), Some("chr2"));

        // Asked for the genotypes alone between two reads, the id of the
        // next variant is neither filled nor left over from this one.
        reader.set_needs(Needs::GTS);
        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.pos, 24);
        assert_eq!(var.gts, [1, 1, 0, MISSING_ALLELE, 0, 0]);
        assert_eq!(var.id, "");
        assert_eq!(var.filled, Needs::GTS | Needs::CHROM_POS);
        assert_eq!(reader.chroms().name(var.chrom), Some("chr1"));

        assert_eq!(reader.chroms().len(), 2);
        assert!(!reader.read_variant(&mut var).unwrap());
        assert!(!reader.read_variant(&mut var).unwrap());
        assert!(var.gts.is_empty());
        assert_eq!(var.filled, Needs::empty());
    }
}
