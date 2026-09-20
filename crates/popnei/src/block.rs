//! Blocks of variants: a run of consecutive variants held as arrays, with
//! the genotypes of all of them in one, and the trait of everything that
//! gives blocks.
//!
//! A block is how the variants flow through popnei, from a source to a
//! calculation: a reader gives blocks, a filter compacts them with
//! [`Block::retain_vars`], a calculation walks their rows with
//! [`Block::variants`] or takes them as a matrix, and a block is also the
//! only way genotypes leave the core, for the Python or the TypeScript user
//! who asks for them. [`BlockReader`] is what every reader implements, and
//! [`Reblock`] is the reader over a reader that cuts and joins the blocks
//! of its source to one size.
//!
//! [`BlockCollector`], which builds blocks by copying the variants that a
//! reader of single variants gives one at a time, and [`CollectedBlocks`],
//! which gives its blocks as a [`BlockReader`], go out with the single
//! variant in work package 3 of `docs/plans/block-readers.md`.
//!
//! `docs/specs/block.md` has the design and section 2 of
//! `docs/architecture.md` the reasons for the arrays.

use std::collections::TryReserveError;
use std::fmt;

use crate::error::{Error, Result};
use crate::variant::{ChromTable, Needs, Variant, VariantReader, VariantRef};

/// How many genotypes a block holds when the caller asks for no number of
/// variants, 5 million, which is `DEF_NUM_GTS_PER_CHUNK` of pyNei's
/// `config.py`. At the ploidy 2 they are the 10 MB of one allocation.
/// Nobody has measured it for popnei.
pub const GENOTYPES_PER_BLOCK: usize = 5_000_000;

/// The fewest variants a block holds when the caller asks for no number,
/// which is what decides for a dataset of many individuals: 100, pyNei's
/// `MIN_NUM_VARS_PER_CHUNK`, measured for popnei by nobody.
pub const MIN_NUM_VARS_PER_BLOCK: usize = 100;

/// The most variants a block holds when the caller asks for no number,
/// which is what decides for a dataset of few individuals: 10000, pyNei's
/// `MAX_NUM_VARS_PER_CHUNK`, measured for popnei by nobody.
pub const MAX_NUM_VARS_PER_BLOCK: usize = 10_000;

/// How many variants a block holds for that many individuals when the
/// caller asks for no number: [`GENOTYPES_PER_BLOCK`] divided by the
/// individuals, never below [`MIN_NUM_VARS_PER_BLOCK`] and never above
/// [`MAX_NUM_VARS_PER_BLOCK`].
///
/// A source of no individual gives [`MAX_NUM_VARS_PER_BLOCK`], which is
/// what the division by one individual would give too.
///
/// It is pyNei's `calc_num_vars_per_chunk` of `pynei/variants.py`: a block
/// is sized by the genotypes it holds and not by its variants, because
/// that is what the memory and the work depend on. pyNei divides by
/// `max(num_samples, 1)`, and so gives the same for no individual.
#[must_use]
pub fn default_num_vars_per_block(num_individuals: usize) -> usize {
    let num_vars = GENOTYPES_PER_BLOCK
        .checked_div(num_individuals)
        .unwrap_or(MAX_NUM_VARS_PER_BLOCK);
    num_vars.clamp(MIN_NUM_VARS_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK)
}

/// The name of each column of a block in Python and in TypeScript, in the
/// order of the columns of [`Block`]. The genotypes are not among them:
/// every block holds them, and no user asks for them.
pub const FIELD_NAMES: [&str; 5] = ["chrom", "pos", "id", "alleles", "qual"];

/// The fields of a variant that each of those names asks the reader for.
/// The chromosome and the position travel together in the core, so either
/// name fills both.
const FIELDS_OF_THE_NAMES: [(&str, Needs); 5] = [
    ("chrom", Needs::CHROM_POS),
    ("pos", Needs::CHROM_POS),
    ("id", Needs::ID),
    ("alleles", Needs::ALLELES),
    ("qual", Needs::QUAL),
];

/// The fields a collector has to be asked for to fill the columns that
/// `names` name, the genotypes among them, which every block holds.
///
/// It is what the binding crates call with the names their user gave, so
/// that one list of names serves Python and TypeScript and a column added
/// later cannot reach one language and not the other.
///
/// # Errors
///
/// When a name is not one of [`FIELD_NAMES`]: the error names it and lists
/// the five, and it is the `ValueError` of Python and the `Error` of
/// TypeScript.
pub fn needs_of_the_fields<'a>(names: impl IntoIterator<Item = &'a str>) -> Result<Needs> {
    let mut needs = Needs::GTS;
    for name in names {
        let found = FIELDS_OF_THE_NAMES
            .iter()
            .find(|(of_a_field, _)| *of_a_field == name);
        let Some((_, field)) = found else {
            return Err(Error::NotAFieldOfABlock {
                name: name.to_string(),
            });
        };
        needs = needs.union(*field);
    }
    Ok(needs)
}

/// The names of the fields of a block between backticks, `` `chrom`,
/// `pos` ``, for the message of a name that is not one of them.
pub(crate) fn field_names_listed() -> String {
    FIELD_NAMES.map(|name| format!("`{name}`")).join(", ")
}

/// How many alleles the buffers of an alleles column are reserved for in
/// each variant, and how many bytes for each of those alleles: the `A` and
/// the `T` of a biallelic SNP, which is the commonest variant of a VCF. A
/// block whose variants have more alleles, or longer ones, grows its
/// buffers from there. Nobody has measured either number.
const ALLELES_PER_VARIANT: usize = 2;
const BYTES_PER_ALLELE: usize = 1;

/// The alleles of the variants of one block, the reference allele of each
/// variant first and then its alternative ones, as the text the source
/// gave.
///
/// The texts are one buffer with the end of each allele in it, and not a
/// string for each allele: a block of 10000 variants is 10000 allocations
/// either way, and one buffer makes it one.
#[derive(Debug)]
pub struct AllelesColumn {
    /// The text of every allele of the block, one after another, as the
    /// bytes of UTF-8: compacting the column moves the bytes of a variant
    /// over the ones of a variant that was dropped, which a `String` does
    /// not do without `unsafe`.
    texts: Vec<u8>,
    /// Where each allele ends in `texts`, one number for each allele.
    allele_ends: Vec<usize>,
    /// Where the alleles of each variant end in `allele_ends`, one number
    /// for each variant.
    var_ends: Vec<usize>,
}

impl AllelesColumn {
    /// A column with no variant in it, whose three buffers are reserved for
    /// `num_vars` variants of [`ALLELES_PER_VARIANT`] alleles of
    /// [`BYTES_PER_ALLELE`] bytes.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the buffers, which the
    /// collector turns into the error of the crate that names the size that
    /// was asked for.
    fn with_num_vars(num_vars: usize) -> std::result::Result<AllelesColumn, TryReserveError> {
        // A number of variants that makes these saturate is far above what
        // any machine gives, and the reservation of it is the error.
        let num_alleles = num_vars.saturating_mul(ALLELES_PER_VARIANT);
        let num_bytes = num_alleles.saturating_mul(BYTES_PER_ALLELE);
        let mut texts = Vec::new();
        texts.try_reserve_exact(num_bytes)?;
        let mut allele_ends = Vec::new();
        allele_ends.try_reserve_exact(num_alleles)?;
        let mut var_ends = Vec::new();
        var_ends.try_reserve_exact(num_vars)?;
        Ok(AllelesColumn {
            texts,
            allele_ends,
            var_ends,
        })
    }

    /// The alleles of one more variant, at the end of the column.
    fn push(&mut self, alleles: &[String]) {
        for allele in alleles {
            self.texts.extend_from_slice(allele.as_bytes());
            self.allele_ends.push(self.texts.len());
        }
        self.var_ends.push(self.allele_ends.len());
    }

    /// Where the alleles of the variants before `var` end in `texts`: the
    /// first byte of the variant `var`, and the length of `texts` for a
    /// `var` beyond the column.
    fn byte_before(&self, var: usize) -> usize {
        let Some(allele) = var
            .checked_sub(1)
            .and_then(|before| self.var_ends.get(before).copied())
        else {
            return 0;
        };
        allele
            .checked_sub(1)
            .and_then(|before| self.allele_ends.get(before).copied())
            .unwrap_or(0)
    }

    /// The variants whose `keep` is true, in their order, and the rest
    /// dropped. The buffers keep their capacity: the bytes and the ends of
    /// a variant that stays move over the ones of a variant that goes.
    fn retain_vars(&mut self, keep: &[bool]) {
        let mut write_byte = 0usize;
        let mut write_allele = 0usize;
        let mut write_var = 0usize;
        let mut read_byte = 0usize;
        let mut read_allele = 0usize;
        for (var, keep_it) in keep.iter().enumerate() {
            let Some(allele_end) = self.var_ends.get(var).copied() else {
                break;
            };
            let byte_end = allele_end
                .checked_sub(1)
                .and_then(|last| self.allele_ends.get(last).copied())
                .unwrap_or(read_byte);
            if *keep_it {
                // How far back the bytes of this variant move, which is
                // what every end of one of its alleles loses. The write
                // position is never past the read one, so the copy never
                // overwrites bytes that have not been moved yet.
                let shift = read_byte.saturating_sub(write_byte);
                if read_byte <= byte_end && byte_end <= self.texts.len() {
                    self.texts.copy_within(read_byte..byte_end, write_byte);
                }
                for index in read_allele..allele_end {
                    let Some(end) = self.allele_ends.get(index).copied() else {
                        break;
                    };
                    if let Some(slot) = self.allele_ends.get_mut(write_allele) {
                        *slot = end.saturating_sub(shift);
                    }
                    write_allele = write_allele.saturating_add(1);
                }
                write_byte = write_byte.saturating_add(byte_end.saturating_sub(read_byte));
                if let Some(slot) = self.var_ends.get_mut(write_var) {
                    *slot = write_allele;
                }
                write_var = write_var.saturating_add(1);
            }
            read_byte = byte_end;
            read_allele = allele_end;
        }
        self.texts.truncate(write_byte);
        self.allele_ends.truncate(write_allele);
        self.var_ends.truncate(write_var);
    }

    /// The alleles of `other` after the ones this column holds, which is
    /// what joining two blocks does with their alleles.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the three buffers.
    fn try_append(&mut self, other: &AllelesColumn) -> std::result::Result<(), TryReserveError> {
        // Two buffers that are both allocated have lengths that add
        // without overflow, so these saturating additions never saturate.
        let byte_base = self.texts.len();
        let allele_base = self.allele_ends.len();
        self.texts.try_reserve(other.texts.len())?;
        self.allele_ends.try_reserve(other.allele_ends.len())?;
        self.var_ends.try_reserve(other.var_ends.len())?;
        self.texts.extend_from_slice(&other.texts);
        self.allele_ends.extend(
            other
                .allele_ends
                .iter()
                .map(|end| end.saturating_add(byte_base)),
        );
        self.var_ends.extend(
            other
                .var_ends
                .iter()
                .map(|end| end.saturating_add(allele_base)),
        );
        Ok(())
    }

    /// The alleles of the variants from `at` on, taken out of this column
    /// into a new one, which is what cutting a block does with them. This
    /// column keeps the first `at` variants and its capacity.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the new column.
    fn try_split_off(&mut self, at: usize) -> std::result::Result<AllelesColumn, TryReserveError> {
        let at = at.min(self.var_ends.len());
        let allele_at = at
            .checked_sub(1)
            .and_then(|before| self.var_ends.get(before).copied())
            .unwrap_or(0);
        let byte_at = self.byte_before(at);
        let mut texts = Vec::new();
        texts.try_reserve_exact(self.texts.len().saturating_sub(byte_at))?;
        texts.extend_from_slice(self.texts.get(byte_at..).unwrap_or(&[]));
        let mut allele_ends = Vec::new();
        allele_ends.try_reserve_exact(self.allele_ends.len().saturating_sub(allele_at))?;
        allele_ends.extend(
            self.allele_ends
                .get(allele_at..)
                .unwrap_or(&[])
                .iter()
                .map(|end| end.saturating_sub(byte_at)),
        );
        let mut var_ends = Vec::new();
        var_ends.try_reserve_exact(self.var_ends.len().saturating_sub(at))?;
        var_ends.extend(
            self.var_ends
                .get(at..)
                .unwrap_or(&[])
                .iter()
                .map(|end| end.saturating_sub(allele_at)),
        );
        self.texts.truncate(byte_at);
        self.allele_ends.truncate(allele_at);
        self.var_ends.truncate(at);
        Ok(AllelesColumn {
            texts,
            allele_ends,
            var_ends,
        })
    }

    /// Where the alleles of `var` are in `allele_ends`, the first and the
    /// one after the last. Both are 0 for a variant that is not there.
    fn alleles_of(&self, var: usize) -> (usize, usize) {
        let Some(end) = self.var_ends.get(var).copied() else {
            return (0, 0);
        };
        let start = var
            .checked_sub(1)
            .and_then(|before| self.var_ends.get(before).copied())
            .unwrap_or(0);
        (start, end)
    }

    /// How many variants the column holds.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.var_ends.len()
    }

    /// How many alleles the variant `var` of the block has, the reference
    /// allele among them. 0 for a variant the column does not hold.
    #[must_use]
    pub fn num_alleles(&self, var: usize) -> usize {
        let (start, end) = self.alleles_of(var);
        end.saturating_sub(start)
    }

    /// The text of the allele `allele` of the variant `var`, as the source
    /// gave it: `A`, `<DEL>`, `*`. The allele 0 is the reference one.
    ///
    /// A variant or an allele the column does not hold gives an empty
    /// text, which no allele of a source is.
    #[must_use]
    pub fn allele(&self, var: usize, allele: usize) -> &str {
        let (first, end) = self.alleles_of(var);
        let Some(index) = first.checked_add(allele).filter(|index| *index < end) else {
            return "";
        };
        let Some(text_end) = self.allele_ends.get(index).copied() else {
            return "";
        };
        let text_start = index
            .checked_sub(1)
            .and_then(|before| self.allele_ends.get(before).copied())
            .unwrap_or(0);
        self.texts
            .get(text_start..text_end)
            .and_then(|bytes| std::str::from_utf8(bytes).ok())
            .unwrap_or("")
    }
}

/// A run of consecutive variants of one source, held as arrays.
///
/// A column other than the genotypes is there only when the collector was
/// asked for it, and `None` when it was not. The number of a chromosome is
/// a number of the table of the reader the block came from, and that table
/// grows while the source is read, so the name of a number is looked up
/// after the block was collected.
#[derive(Debug)]
pub struct Block {
    /// How many variants the block holds, which is the number that was
    /// asked for except in the last block of a source.
    pub num_vars: usize,
    /// How many individuals the source has, the same for every variant.
    pub num_individuals: usize,
    /// How many alleles the genotype of one individual holds.
    pub ploidy: usize,
    /// `num_vars` x `num_individuals` x `ploidy` alleles, variant after
    /// variant and inside a variant individual after individual. 0 is the
    /// reference allele, 1 and above the alternative ones, and
    /// [`MISSING_ALLELE`](crate::variant::MISSING_ALLELE), -1, an allele
    /// that was not called.
    pub gts: Vec<i8>,
    /// The number of the chromosome of each variant, in the
    /// [`ChromTable`](crate::variant::ChromTable) of the reader the block
    /// came from.
    pub chrom: Option<Vec<u32>>,
    /// The position of each variant, 1 based as in a VCF.
    pub pos: Option<Vec<u64>>,
    /// The id of each variant, empty for a variant that has none.
    pub id: Option<Vec<String>>,
    /// The alleles of each variant, the reference one first.
    pub alleles: Option<AllelesColumn>,
    /// The quality of each variant, phred scaled as the QUAL of a VCF, and
    /// NaN for a variant that has none.
    pub qual: Option<Vec<f32>>,
}

impl Block {
    /// Which fields the block holds: the columns that are there, and the
    /// genotypes when `gts` is not empty or the block has no variant.
    ///
    /// A consumer that depends on a field asks the block for its fields and
    /// fails with the error that names what is missing,
    /// `asked_for.difference(block.fields())`.
    #[must_use]
    pub fn fields(&self) -> Needs {
        let mut fields = Needs::empty();
        if !self.gts.is_empty() || self.num_vars == 0 {
            fields |= Needs::GTS;
        }
        // The chromosome and the position are one field, and a block that
        // holds one of the two columns and not the other is a block that
        // `check` does not look at: it holds neither field.
        if self.chrom.is_some() && self.pos.is_some() {
            fields |= Needs::CHROM_POS;
        }
        if self.id.is_some() {
            fields |= Needs::ID;
        }
        if self.alleles.is_some() {
            fields |= Needs::ALLELES;
        }
        if self.qual.is_some() {
            fields |= Needs::QUAL;
        }
        fields
    }

    /// Which columns the block has, one flag for the genotypes and one for
    /// each of the five, so that two blocks are joined only when every
    /// column of the one is a column of the other.
    ///
    /// It is not [`Block::fields`]: that one answers what a consumer can
    /// read, and puts the chromosome and the position together.
    fn columns(&self) -> [bool; 6] {
        [
            !self.gts.is_empty(),
            self.chrom.is_some(),
            self.pos.is_some(),
            self.id.is_some(),
            self.alleles.is_some(),
            self.qual.is_some(),
        ]
    }

    /// How many alleles one variant of the block holds, the individuals
    /// times the ploidy, or the error of a block that cannot be.
    fn gts_per_var(&self) -> Result<usize> {
        self.num_individuals
            .checked_mul(self.ploidy)
            .ok_or(Error::BlockTooLarge {
                num_vars_per_block: self.num_vars,
                num_individuals: self.num_individuals,
                ploidy: self.ploidy,
            })
    }

    /// The views of the variants of the block, in order.
    ///
    /// A view allocates nothing: its genotypes are a slice of `gts` and its
    /// other fields are read out of the columns. A block whose arrays are
    /// not of its size, which [`Block::check`] finds, gives the views up to
    /// the first variant that is not in them.
    pub fn variants(&self) -> impl Iterator<Item = VariantRef<'_>> {
        (0..self.num_vars).map_while(|var| self.variant(var))
    }

    /// The view of the variant `var` of the block, counted from 0, or
    /// `None` when the block has no such variant.
    #[must_use]
    pub fn variant(&self, var: usize) -> Option<VariantRef<'_>> {
        if var >= self.num_vars {
            return None;
        }
        let gts = match self.gts.is_empty() {
            // A block built without the genotypes: every view has none.
            true => &[][..],
            false => {
                let gts_per_var = self.gts_per_var().ok()?;
                let start = var.checked_mul(gts_per_var)?;
                let end = start.checked_add(gts_per_var)?;
                self.gts.get(start..end)?
            }
        };
        Some(VariantRef::new(
            gts,
            self.chrom
                .as_ref()
                .and_then(|column| column.get(var))
                .copied(),
            self.pos
                .as_ref()
                .and_then(|column| column.get(var))
                .copied(),
            self.id
                .as_ref()
                .and_then(|column| column.get(var))
                .map(String::as_str),
            self.qual
                .as_ref()
                .and_then(|column| column.get(var))
                .copied(),
            self.alleles.as_ref().map(|column| (column, var)),
        ))
    }

    /// It keeps the variants whose `keep` is true, in their order, in the
    /// genotypes and in every column, in place, and sets `num_vars` to how
    /// many stayed. It is what a filter of variants calls with the rows it
    /// decided to keep.
    ///
    /// Nothing is allocated and the block keeps its capacity. A `gts` that
    /// is empty, of a block built without the genotypes, stays empty.
    ///
    /// # Errors
    ///
    /// When `keep` has not one value for each variant of the block, and
    /// when the arrays of the block are not of its size, which
    /// [`Block::check`] finds. In both the block is left as it was.
    pub fn retain_vars(&mut self, keep: &[bool]) -> Result<()> {
        if keep.len() != self.num_vars {
            return Err(Error::KeepOfAnotherSize {
                found: keep.len(),
                num_vars: self.num_vars,
            });
        }
        // The rows are moved by their place in `gts`, so the arrays have to
        // be of the size the block says before any of them is touched.
        self.check()?;
        let gts_per_var = self.gts_per_var()?;
        if !self.gts.is_empty() {
            let mut write = 0usize;
            for (var, keep_it) in keep.iter().enumerate() {
                if !*keep_it {
                    continue;
                }
                // Every one of these is a place in `gts`, which `check`
                // just said holds num_vars x gts_per_var alleles.
                let read = var.saturating_mul(gts_per_var);
                let end = read.saturating_add(gts_per_var);
                if end <= self.gts.len() {
                    self.gts.copy_within(read..end, write);
                }
                write = write.saturating_add(gts_per_var);
            }
            self.gts.truncate(write);
        }
        retain_in_column(self.chrom.as_mut(), keep);
        retain_in_column(self.pos.as_mut(), keep);
        retain_in_column(self.id.as_mut(), keep);
        retain_in_column(self.qual.as_mut(), keep);
        if let Some(alleles) = self.alleles.as_mut() {
            alleles.retain_vars(keep);
        }
        self.num_vars = keep.iter().filter(|keep_it| **keep_it).count();
        Ok(())
    }

    /// That `gts` holds `num_vars` x `num_individuals` x `ploidy` alleles,
    /// or none, and that every column that is there holds `num_vars`
    /// entries.
    ///
    /// The fields of a block are public, so a reader with a defect can
    /// build one that breaks this, and its genotypes would then be read
    /// one at the place of another with nothing to show it. It is called
    /// where that would happen: by the binding crates before the genotypes
    /// cross to numpy or to an `Int8Array`, by [`Reblock`] on every block
    /// it takes, and by the writer of a vars file.
    ///
    /// # Errors
    ///
    /// When an array is not of the size of the block: the error names the
    /// array, how many entries it holds and how many the block says.
    pub fn check(&self) -> Result<()> {
        let gts_per_var = self.gts_per_var()?;
        let gts_of_the_block =
            self.num_vars
                .checked_mul(gts_per_var)
                .ok_or(Error::BlockTooLarge {
                    num_vars_per_block: self.num_vars,
                    num_individuals: self.num_individuals,
                    ploidy: self.ploidy,
                })?;
        if !self.gts.is_empty() && self.gts.len() != gts_of_the_block {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "gts",
                found: self.gts.len(),
                expected: gts_of_the_block,
            });
        }
        let lengths = [
            ("chrom", self.chrom.as_ref().map(Vec::len)),
            ("pos", self.pos.as_ref().map(Vec::len)),
            ("id", self.id.as_ref().map(Vec::len)),
            ("qual", self.qual.as_ref().map(Vec::len)),
            (
                "alleles",
                self.alleles.as_ref().map(AllelesColumn::num_vars),
            ),
        ];
        for (array, length) in lengths {
            if let Some(length) = length
                && length != self.num_vars
            {
                return Err(Error::BlockArrayOfAnotherSize {
                    array,
                    found: length,
                    expected: self.num_vars,
                });
            }
        }
        Ok(())
    }
}

/// The entries of one column whose `keep` is true, in their order, and the
/// rest dropped. A column that is not in the block is left alone.
///
/// `Vec::retain` visits the entries in order and moves the ones that stay
/// down over the ones that go, so nothing is allocated.
fn retain_in_column<T>(column: Option<&mut Vec<T>>, keep: &[bool]) {
    let Some(column) = column else {
        return;
    };
    let mut keep = keep.iter();
    column.retain(|_| keep.next().copied().unwrap_or(false));
}

/// Anything that gives blocks: the VCF reader, the vars file reader, a
/// filter over another reader, [`Reblock`].
///
/// The contract, which each reader keeps itself:
///
/// - A block it gives holds one variant at least, and the block is the
///   caller's.
/// - `next_block` gives `None` when the source has no more variants, and
///   `None` at every call after that.
/// - After an error it gives `None` at every call, and a reader over
///   another reader does not call its source again. A reader that went on
///   would give the variants that follow a wrong one as if nothing had
///   happened.
///
/// It can be used as a boxed trait object, `Box<dyn BlockReader>`, which is
/// how the two binding crates hold their reader, because neither a pyo3
/// class nor a wasm-bindgen class can be generic; so it has no generic
/// method and no method that takes or gives `Self`, and it is implemented
/// for `Box<dyn BlockReader>` too, so that what is generic over a reader,
/// a filter or [`Reblock`], takes a boxed one. It asks for `Send` because
/// a read ahead thread moves a reader into another thread.
pub trait BlockReader: Send {
    /// The next block, which holds one variant at least, or `None` when
    /// there are no more.
    ///
    /// # Errors
    ///
    /// When the source cannot be read or what it holds is malformed. The
    /// block the error happened in is lost, the blocks before it were
    /// given, and every call after it gives `None`.
    fn next_block(&mut self) -> Result<Option<Block>>;

    /// The names of the individuals, in the order of their genotypes in
    /// the rows of a block.
    fn individuals(&self) -> &[String];

    /// How many alleles the genotype of one individual holds.
    fn ploidy(&self) -> usize;

    /// The names of the chromosomes seen so far, each with its number. It
    /// grows while the source is read, so the name of the number of a
    /// variant is looked up after the block was given.
    fn chroms(&self) -> &ChromTable;

    /// Which fields the reader is asked to fill. The rest may be skipped,
    /// and a column that was not asked for is `None` in the block.
    /// [`Needs::ALL`] until it is called, and the change holds from the
    /// next block that is built.
    fn set_needs(&mut self, needs: Needs);
}

impl<R: BlockReader + ?Sized> BlockReader for Box<R> {
    fn next_block(&mut self) -> Result<Option<Block>> {
        (**self).next_block()
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

/// A reader over a reader that gives the variants of its source in blocks
/// of one size, the last one aside: it joins the blocks that are too short
/// and cuts the ones that are too long.
///
/// It goes where the size matters: before the matrix work when a filter
/// took variants out, before the writer of a vars file, and at the end of
/// `iter_blocks`, so that a user gets blocks of the size they asked for. A
/// block of its source that already has the size, with nothing waiting from
/// the one before, goes through as it is, with no copy.
///
/// It keeps at most one block from one call to the next, so its memory is
/// two blocks. Joining copies the rows of the block that arrives after the
/// ones that were waiting, and cutting copies the rows after the cut: one
/// copy of each column per block and none per variant.
pub struct Reblock<R: BlockReader> {
    reader: R,
    num_vars_per_block: usize,
    /// The individuals and the ploidy of the source, which every block it
    /// gives has to have for its rows to be joined with the others'.
    num_individuals: usize,
    ploidy: usize,
    /// The variants that did not fill a block, or the ones left after a
    /// cut. It holds one variant at least when it is there.
    waiting: Option<Block>,
    /// Whether the source has no more blocks or gave an error. After
    /// either there is no block.
    finished: bool,
}

impl<R: BlockReader> Reblock<R> {
    /// The reader that gives the variants of `reader` in blocks of
    /// `num_vars_per_block` variants.
    ///
    /// `num_vars_per_block` is 1 or more, or `None` for
    /// [`default_num_vars_per_block`] for the individuals of `reader`.
    ///
    /// # Errors
    ///
    /// When `num_vars_per_block` is 0, and when the genotypes of one block,
    /// the variants times the individuals times the ploidy, are more than
    /// this machine addresses.
    pub fn new(reader: R, num_vars_per_block: Option<usize>) -> Result<Reblock<R>> {
        let num_individuals = reader.individuals().len();
        let ploidy = reader.ploidy();
        let num_vars_per_block = match num_vars_per_block {
            Some(0) => return Err(Error::BlockOfNoVariants),
            Some(asked_for) => asked_for,
            None => default_num_vars_per_block(num_individuals),
        };
        let too_large = || Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
        };
        let gts_per_var = num_individuals.checked_mul(ploidy).ok_or_else(too_large)?;
        gts_per_var
            .checked_mul(num_vars_per_block)
            .ok_or_else(too_large)?;
        Ok(Reblock {
            reader,
            num_vars_per_block,
            num_individuals,
            ploidy,
            waiting: None,
            finished: false,
        })
    }

    /// The error of a block whose memory the machine does not give.
    fn too_large(&self) -> Error {
        Error::BlockTooLarge {
            num_vars_per_block: self.num_vars_per_block,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
        }
    }

    /// That a block of the source is one whose rows can be joined with the
    /// rows of the blocks before it: its arrays are of its size, it holds a
    /// variant, and its individuals and its ploidy are the ones the source
    /// says it has.
    fn taken(&self, block: &Block) -> Result<()> {
        block.check()?;
        // A source that gives a block of no variants has a defect, and it
        // is not asked again: over a source that always gives one, a
        // `reblock` that asked again would never come back.
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        if block.num_individuals != self.num_individuals || block.ploidy != self.ploidy {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals: self.num_individuals,
                ploidy: self.ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        Ok(())
    }

    /// The rows of `block` after the ones of `waiting`, which have the same
    /// columns: one copy of each column and none per variant.
    fn join(&self, waiting: &mut Block, block: Block) -> Result<()> {
        let num_vars = waiting
            .num_vars
            .checked_add(block.num_vars)
            .ok_or_else(|| self.too_large())?;
        try_extend(&mut waiting.gts, block.gts).map_err(|_| self.too_large())?;
        try_extend_column(waiting.chrom.as_mut(), block.chrom).map_err(|_| self.too_large())?;
        try_extend_column(waiting.pos.as_mut(), block.pos).map_err(|_| self.too_large())?;
        try_extend_column(waiting.id.as_mut(), block.id).map_err(|_| self.too_large())?;
        try_extend_column(waiting.qual.as_mut(), block.qual).map_err(|_| self.too_large())?;
        if let (Some(waiting), Some(block)) = (waiting.alleles.as_mut(), block.alleles.as_ref()) {
            waiting.try_append(block).map_err(|_| self.too_large())?;
        }
        waiting.num_vars = num_vars;
        Ok(())
    }

    /// The block that is waiting, cut to the size that was asked for: the
    /// first `num_vars_per_block` variants are given and the rest wait. A
    /// block that has the size already is given as it is, with no copy.
    fn cut(&mut self, mut waiting: Block) -> Result<Block> {
        if waiting.num_vars <= self.num_vars_per_block {
            return Ok(waiting);
        }
        let left_over = split_block(&mut waiting, self.num_vars_per_block)
            .inspect_err(|_| {
                self.finished = true;
            })
            .map_err(|_| self.too_large())?;
        self.waiting = Some(left_over);
        Ok(waiting)
    }
}

impl<R: BlockReader> BlockReader for Reblock<R> {
    /// The next block of `num_vars_per_block` variants, the last one of the
    /// source aside, which is the only one that can be shorter; or the
    /// shorter block that was waiting when the columns of the source
    /// changed.
    ///
    /// # Errors
    ///
    /// When the source fails, which loses the variants that were waiting
    /// too; when a block of the source is not of its own size, holds no
    /// variant, or has other individuals or another ploidy than the source
    /// says; and when the machine does not give the memory of a block.
    /// After any of them there is no block.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            if let Some(waiting) = self.waiting.take() {
                if waiting.num_vars >= self.num_vars_per_block {
                    return self.cut(waiting).map(Some);
                }
                self.waiting = Some(waiting);
            }
            let block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    // What was waiting is the last block, shorter than the
                    // size that was asked for.
                    return Ok(self.waiting.take());
                }
                Err(error) => {
                    // The source ends at its error and is not called
                    // again, and what was waiting is lost with it.
                    self.finished = true;
                    self.waiting = None;
                    return Err(error);
                }
            };
            if let Err(error) = self.taken(&block) {
                self.finished = true;
                self.waiting = None;
                return Err(error);
            }
            match self.waiting.take() {
                None => self.waiting = Some(block),
                Some(mut waiting) => {
                    if waiting.columns() != block.columns() {
                        // The columns of the source changed, which a
                        // change of `Needs` in the middle of a pass does:
                        // what was waiting is given as a shorter block.
                        self.waiting = Some(block);
                        return Ok(Some(waiting));
                    }
                    if let Err(error) = self.join(&mut waiting, block) {
                        self.finished = true;
                        return Err(error);
                    }
                    self.waiting = Some(waiting);
                }
            }
        }
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    /// The table of the source: a reader over another reader has none of
    /// its own.
    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    /// What the source is asked for. The block that is waiting keeps the
    /// columns it was built with, and it is given on its own, since blocks
    /// of different columns are not joined.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs);
    }
}

impl<R: BlockReader> fmt::Debug for Reblock<R> {
    /// The size it gives and where it has got to. The reader is left out,
    /// so that a `Reblock` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Reblock")
            .field("num_vars_per_block", &self.num_vars_per_block)
            .field("num_individuals", &self.num_individuals)
            .field("ploidy", &self.ploidy)
            .field(
                "waiting",
                &self.waiting.as_ref().map(|block| block.num_vars),
            )
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

/// `values` after what `column` holds, with the memory asked of the machine
/// first. A column that is in neither block, or in one of them alone, is
/// left as it is: `Reblock` joins blocks of the same columns.
fn try_extend_column<T>(
    column: Option<&mut Vec<T>>,
    values: Option<Vec<T>>,
) -> std::result::Result<(), TryReserveError> {
    let (Some(column), Some(values)) = (column, values) else {
        return Ok(());
    };
    try_extend(column, values)
}

/// `values` after what `column` holds, with the memory asked of the machine
/// with `try_reserve`: `Vec::extend` ends the process when the machine has
/// not the memory, and a block of a size that a caller wrote reaches it.
fn try_extend<T>(column: &mut Vec<T>, values: Vec<T>) -> std::result::Result<(), TryReserveError> {
    column.try_reserve(values.len())?;
    column.extend(values);
    Ok(())
}

/// The variants of `block` from `at` on, taken out of it into a new block,
/// which keeps the columns and the individuals of the one it came from.
/// `block` keeps its first `at` variants and its capacity.
///
/// # Errors
///
/// When the machine does not give the memory of the new block.
fn split_block(block: &mut Block, at: usize) -> std::result::Result<Block, TryReserveError> {
    let num_vars = block.num_vars.saturating_sub(at);
    // The block holds `at` + `num_vars` variants, so these places are in
    // its arrays: `check` said so before it was taken.
    let gts_per_var = block.num_individuals.saturating_mul(block.ploidy);
    let mut gts = Vec::new();
    let gts_at = at.saturating_mul(gts_per_var);
    if !block.gts.is_empty() {
        gts.try_reserve_exact(block.gts.len().saturating_sub(gts_at))?;
        gts.extend_from_slice(block.gts.get(gts_at..).unwrap_or(&[]));
        block.gts.truncate(gts_at);
    }
    let alleles = match block.alleles.as_mut() {
        Some(alleles) => Some(alleles.try_split_off(at)?),
        None => None,
    };
    let left_over = Block {
        num_vars,
        num_individuals: block.num_individuals,
        ploidy: block.ploidy,
        gts,
        chrom: try_split_column(block.chrom.as_mut(), at)?,
        pos: try_split_column(block.pos.as_mut(), at)?,
        id: try_split_column(block.id.as_mut(), at)?,
        alleles,
        qual: try_split_column(block.qual.as_mut(), at)?,
    };
    block.num_vars = at;
    Ok(left_over)
}

/// The entries of `column` from `at` on, taken out of it into a new column
/// whose memory is asked for with `try_reserve`. `Vec::split_off` allocates
/// without asking, and a size that a caller wrote reaches it.
fn try_split_column<T>(
    column: Option<&mut Vec<T>>,
    at: usize,
) -> std::result::Result<Option<Vec<T>>, TryReserveError> {
    let Some(column) = column else {
        return Ok(None);
    };
    let at = at.min(column.len());
    let mut left_over = Vec::new();
    left_over.try_reserve_exact(column.len().saturating_sub(at))?;
    left_over.extend(column.drain(at..));
    Ok(Some(left_over))
}

/// It builds blocks from a reader, one after another, each of the number of
/// variants that was asked for.
///
/// It owns its reader and gives it back to whoever wants the individuals,
/// the ploidy or the table of the chromosomes. What it keeps from one block
/// to the next is the [`Variant`] it lends to the reader; each block is
/// allocated when it is started and given away when it is full.
pub struct BlockCollector<R: VariantReader> {
    reader: R,
    /// What each block holds, the genotypes among them.
    needs: Needs,
    num_vars_per_block: usize,
    num_individuals: usize,
    ploidy: usize,
    /// `num_individuals` x `ploidy`, the alleles of one variant, which the
    /// trait asks of every reader and the collector checks in each one.
    gts_per_variant: usize,
    /// `num_vars_per_block` x `num_individuals` x `ploidy`, the genotypes
    /// a full block holds, which every block is allocated for.
    gts_per_block: usize,
    /// The variant the reader is lent, again and again.
    var: Variant,
    /// Whether the reader has no more variants or gave an error. After
    /// either there is no block.
    finished: bool,
}

impl<R: VariantReader> BlockCollector<R> {
    /// The collector over `reader`, which it asks for `needs` and the
    /// genotypes.
    ///
    /// `needs` is what each block will hold besides the genotypes, which
    /// are always part of it, and a field that is not in it is a column
    /// that the block does not have and that the reader is never asked to
    /// parse. `num_vars_per_block` is how many variants a block holds, 1 or
    /// more, and `None` is [`default_num_vars_per_block`] for the
    /// individuals of the reader.
    ///
    /// # Errors
    ///
    /// When `num_vars_per_block` is 0, and when the genotypes of one block,
    /// the variants times the individuals times the ploidy, are more than
    /// the machine addresses.
    pub fn new(mut reader: R, needs: Needs, num_vars_per_block: Option<usize>) -> Result<Self> {
        let num_individuals = reader.individuals().len();
        let ploidy = reader.ploidy();
        let num_vars_per_block = match num_vars_per_block {
            Some(0) => return Err(Error::BlockOfNoVariants),
            Some(asked_for) => asked_for,
            None => default_num_vars_per_block(num_individuals),
        };
        let too_large = || Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
        };
        let gts_per_variant = num_individuals.checked_mul(ploidy).ok_or_else(too_large)?;
        let gts_per_block = gts_per_variant
            .checked_mul(num_vars_per_block)
            .ok_or_else(too_large)?;
        let needs = needs.union(Needs::GTS);
        reader.set_needs(needs);
        Ok(BlockCollector {
            reader,
            needs,
            num_vars_per_block,
            num_individuals,
            ploidy,
            gts_per_variant,
            gts_per_block,
            var: Variant::new(),
            finished: false,
        })
    }

    /// The next block, or `None` when the reader has no more variants.
    ///
    /// The last block of a source is the only one that can hold fewer
    /// variants than were asked for, and a reader with no variants gives no
    /// block at all.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the columns of the
    /// block, which is asked for before a variant is read; when the reader
    /// fails; when it did not fill a field that was asked for; and when it
    /// filled a variant with a number of alleles other than its individuals
    /// times its ploidy. The block that was being built is lost with the
    /// error, and every call after it gives no block.
    pub fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        // The memory of the columns is asked for before the first variant
        // is read, so that a size the caller wrote is refused before the
        // source is touched. The count of the variants is the one of the
        // range, so that a block of `usize::MAX` variants needs no addition
        // of our own.
        let mut block = self.start_block().inspect_err(|_| {
            self.finished = true;
        })?;
        for count in 1..=self.num_vars_per_block {
            // The trait says that a reader ends at its error, and the
            // collector ends too: what a reader that goes on gives after
            // one is not the source, and no caller of ours reads it.
            let read = self.reader.read_variant(&mut self.var).inspect_err(|_| {
                self.finished = true;
            })?;
            if !read {
                self.finished = true;
                break;
            }
            let not_filled = self.needs.difference(self.var.filled);
            if !not_filled.is_empty() {
                self.finished = true;
                return Err(Error::FieldsNotFilled { fields: not_filled });
            }
            // The genotypes of a block are read as variants x individuals
            // x ploidy, and a variant of another number of alleles would
            // move every genotype after it without a sign of it.
            if self.var.gts.len() != self.gts_per_variant {
                self.finished = true;
                return Err(Error::VariantOfAnotherSize {
                    found: self.var.gts.len(),
                    expected: self.gts_per_variant,
                    num_individuals: self.num_individuals,
                    ploidy: self.ploidy,
                });
            }
            self.push_variant(&mut block);
            block.num_vars = count;
        }
        if block.num_vars == 0 {
            return Ok(None);
        }
        Ok(Some(block))
    }

    /// The reader the collector was built over, for the individuals, the
    /// ploidy and the names of the chromosome numbers of a block.
    pub fn reader(&self) -> &R {
        &self.reader
    }

    /// What the blocks from the next one on will hold, and what the reader
    /// is asked for. The genotypes are part of it whatever is given: a
    /// collector copies the genotypes of every variant it reads.
    pub fn set_needs(&mut self, needs: Needs) {
        self.needs = needs.union(Needs::GTS);
        self.reader.set_needs(self.needs);
    }

    /// The error of a block whose memory the machine does not give, with
    /// the size that was asked for.
    fn too_large(&self) -> Error {
        Error::BlockTooLarge {
            num_vars_per_block: self.num_vars_per_block,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
        }
    }

    /// An empty column reserved for `num_items`, or `None` when `field` is
    /// not among what the blocks of this collector hold.
    ///
    /// The memory is asked for with `try_reserve_exact`, which gives it
    /// back as an error: `Vec::with_capacity` ends the process when the
    /// machine has not the memory, and panics above what a `Vec` holds,
    /// and a size that a caller of popnei wrote reaches both.
    fn reserved_column<T>(&self, field: Needs, num_items: usize) -> Result<Option<Vec<T>>> {
        if !self.needs.contains(field) {
            return Ok(None);
        }
        let mut column = Vec::new();
        column
            .try_reserve_exact(num_items)
            .map_err(|_| self.too_large())?;
        Ok(Some(column))
    }

    /// An empty block with every column the collector was asked for, each
    /// reserved for a full block.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of one of the columns.
    fn start_block(&self) -> Result<Block> {
        let num_vars = self.num_vars_per_block;
        let mut gts = Vec::new();
        gts.try_reserve_exact(self.gts_per_block)
            .map_err(|_| self.too_large())?;
        let alleles = match self.needs.contains(Needs::ALLELES) {
            true => Some(AllelesColumn::with_num_vars(num_vars).map_err(|_| self.too_large())?),
            false => None,
        };
        Ok(Block {
            num_vars: 0,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
            gts,
            chrom: self.reserved_column(Needs::CHROM_POS, num_vars)?,
            pos: self.reserved_column(Needs::CHROM_POS, num_vars)?,
            id: self.reserved_column(Needs::ID, num_vars)?,
            alleles,
            qual: self.reserved_column(Needs::QUAL, num_vars)?,
        })
    }

    /// The variant the reader last filled, copied into the columns of the
    /// block. The reader gives the `num_individuals` x `ploidy` alleles of
    /// every variant that the trait asks of it, so the genotypes are one
    /// copy of a few kilobytes.
    fn push_variant(&self, block: &mut Block) {
        block.gts.extend_from_slice(&self.var.gts);
        if let Some(chrom) = block.chrom.as_mut() {
            chrom.push(self.var.chrom);
        }
        if let Some(pos) = block.pos.as_mut() {
            pos.push(self.var.pos);
        }
        if let Some(id) = block.id.as_mut() {
            id.push(self.var.id.clone());
        }
        if let Some(alleles) = block.alleles.as_mut() {
            alleles.push(&self.var.alleles);
        }
        if let Some(qual) = block.qual.as_mut() {
            // A variant with no quality is a NaN in the column, which is
            // what Python and TypeScript are given for it.
            qual.push(self.var.qual.unwrap_or(f32::NAN));
        }
    }
}

/// The blocks a [`BlockCollector`] builds, as a [`BlockReader`].
///
/// It goes out in work package 3 of `docs/plans/block-readers.md`, with the
/// single variant and the collector: it is here so that everything that
/// takes a reader of blocks, the binding crates among them, can be written
/// and tested before the VCF reader gives blocks itself.
///
/// What it does not keep of the contract of [`BlockReader`]: a block of it
/// always holds the genotypes, because a collector copies the genotypes of
/// every variant it reads, so a `set_needs` without
/// [`Needs::GTS`](crate::variant::Needs::GTS) still gives them.
pub struct CollectedBlocks<R: VariantReader> {
    collector: BlockCollector<R>,
}

impl<R: VariantReader> CollectedBlocks<R> {
    /// The blocks of `reader`, each of `num_vars_per_block` variants and
    /// holding `needs` and the genotypes.
    ///
    /// # Errors
    ///
    /// The ones of [`BlockCollector::new`]: a `num_vars_per_block` of 0,
    /// and a block of more genotypes than this machine addresses.
    pub fn new(reader: R, needs: Needs, num_vars_per_block: Option<usize>) -> Result<Self> {
        Ok(CollectedBlocks {
            collector: BlockCollector::new(reader, needs, num_vars_per_block)?,
        })
    }
}

impl<R: VariantReader> BlockReader for CollectedBlocks<R> {
    /// The next block of the collector, which holds one variant at least,
    /// gives `None` when the reader has no more and `None` at every call
    /// after an error.
    ///
    /// # Errors
    ///
    /// The ones of [`BlockCollector::next_block`]: the reader failed, it
    /// did not fill a field that was asked for, it filled a variant of
    /// another number of alleles, or the machine did not give the memory
    /// of a column.
    fn next_block(&mut self) -> Result<Option<Block>> {
        self.collector.next_block()
    }

    fn individuals(&self) -> &[String] {
        self.collector.reader().individuals()
    }

    fn ploidy(&self) -> usize {
        self.collector.reader().ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.collector.reader().chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.collector.set_needs(needs);
    }
}

impl<R: VariantReader> fmt::Debug for CollectedBlocks<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CollectedBlocks")
            .field("collector", &self.collector)
            .finish()
    }
}

impl<R: VariantReader> fmt::Debug for BlockCollector<R> {
    /// What the collector was asked for and where it has got to. The reader
    /// is left out, so that a collector over a reader that has no `Debug`
    /// has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("BlockCollector")
            .field("needs", &self.needs)
            .field("num_vars_per_block", &self.num_vars_per_block)
            .field("num_individuals", &self.num_individuals)
            .field("ploidy", &self.ploidy)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{
        AllelesColumn, Block, BlockCollector, BlockReader, CollectedBlocks, FIELD_NAMES,
        FIELDS_OF_THE_NAMES, GENOTYPES_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK, MIN_NUM_VARS_PER_BLOCK,
        Reblock, default_num_vars_per_block, needs_of_the_fields,
    };
    use crate::error::{Error, Result};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs, Variant, VariantReader, VariantRef};

    /// The reference VCFs live at the root of the repository, beside the
    /// Python tests that read the same files, and not inside this crate.
    /// The path is built from the directory of the manifest, so it holds
    /// whether the tests are run with `cargo test --workspace` or with
    /// `cargo test -p popnei`.
    fn reference_vcf(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/vcf")
            .join(name)
    }

    /// The 50 individuals of `many.vcf`, which
    /// `tests/reference/vcf/make_reference.py` writes as diploid.
    const MANY_INDIVIDUALS: usize = 50;
    const MANY_PLOIDY: usize = 2;

    /// The options of a VCF read with every variant given, the ones that
    /// failed a filter among them.
    fn every_variant() -> VcfOptions {
        VcfOptions {
            ploidy: MANY_PLOIDY,
            only_passed: false,
        }
    }

    /// A collector over one of the reference VCFs.
    fn collector_over(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> BlockCollector<VcfReader<BufReader<File>>> {
        let reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        match BlockCollector::new(reader, needs, num_vars_per_block) {
            Ok(collector) => collector,
            Err(error) => panic!("{name}: the collector was not built: {error}"),
        }
    }

    /// Every block a collector gives, until it has no more or it fails.
    fn blocks_of<R: VariantReader>(collector: &mut BlockCollector<R>) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = collector.next_block()? {
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// The blocks of one of the reference VCFs.
    fn blocks_read(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> Vec<Block> {
        let mut collector = collector_over(name, options, needs, num_vars_per_block);
        match blocks_of(&mut collector) {
            Ok(blocks) => blocks,
            Err(error) => panic!("{name}: the collector stopped at {error}"),
        }
    }

    /// How many variants each block holds.
    fn num_vars_of(blocks: &[Block]) -> Vec<usize> {
        blocks.iter().map(|block| block.num_vars).collect()
    }

    /// How many genotypes each block holds, which the consumer that
    /// reshapes `gts` into variants x individuals x ploidy depends on.
    fn num_gts_of(blocks: &[Block]) -> Vec<usize> {
        blocks.iter().map(|block| block.gts.len()).collect()
    }

    #[test]
    fn the_blocks_of_a_hundred_variants_of_many_vcf_are_four_full_and_one_of_seventy_five() {
        // The default of the reader gives 475 of the 500 variants.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(100));
        assert_eq!(num_vars_of(&blocks), [100, 100, 100, 100, 75]);
        // 100 x 50 x 2 genotypes in a full block and 75 x 50 x 2 in the
        // last one.
        assert_eq!(num_gts_of(&blocks), [10_000, 10_000, 10_000, 10_000, 7_500]);
        for block in &blocks {
            assert_eq!(block.num_individuals, MANY_INDIVIDUALS);
            assert_eq!(block.ploidy, MANY_PLOIDY);
        }
    }

    #[test]
    fn every_variant_of_many_vcf_given_makes_five_blocks_of_a_hundred() {
        let blocks = blocks_read("many.vcf", every_variant(), Needs::GTS, Some(100));
        assert_eq!(num_vars_of(&blocks), [100, 100, 100, 100, 100]);
        assert_eq!(num_gts_of(&blocks), [10_000; 5]);
    }

    #[test]
    fn a_block_of_a_thousand_variants_holds_the_four_hundred_and_seventy_five_of_the_default() {
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(1000));
        assert_eq!(num_vars_of(&blocks), [475]);
        // 475 x 50 x 2.
        assert_eq!(num_gts_of(&blocks), [47_500]);
    }

    /// One variant as `read_variant` gives it, which is what the blocks,
    /// joined, have to hold whatever their size.
    #[derive(Debug, PartialEq, Eq)]
    struct Site {
        chrom: u32,
        pos: u64,
        gts: Vec<i8>,
    }

    /// The variants of a VCF read one by one, which is what the blocks are
    /// compared with.
    fn sites_read_one_by_one(name: &str, options: VcfOptions) -> Vec<Site> {
        let mut reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let mut var = Variant::new();
        let mut sites = Vec::new();
        loop {
            match reader.read_variant(&mut var) {
                Ok(true) => sites.push(Site {
                    chrom: var.chrom,
                    pos: var.pos,
                    gts: var.gts.clone(),
                }),
                Ok(false) => return sites,
                Err(error) => panic!("{name}: the reader stopped at {error}"),
            }
        }
    }

    /// The variants of the blocks, joined, one by one.
    fn sites_of(blocks: &[Block]) -> Vec<Site> {
        let mut sites = Vec::new();
        for block in blocks {
            let chrom = block
                .chrom
                .as_ref()
                .expect("the chromosomes were asked for");
            let pos = block.pos.as_ref().expect("the positions were asked for");
            let gts_per_var = block.num_individuals.saturating_mul(block.ploidy);
            assert_eq!(chrom.len(), block.num_vars);
            assert_eq!(pos.len(), block.num_vars);
            assert_eq!(block.gts.len(), block.num_vars.saturating_mul(gts_per_var));
            for (index, gts) in block.gts.chunks_exact(gts_per_var).enumerate() {
                sites.push(Site {
                    chrom: chrom[index],
                    pos: pos[index],
                    gts: gts.to_vec(),
                });
            }
        }
        sites
    }

    #[test]
    fn the_blocks_joined_are_the_variants_the_reader_gives_one_by_one() {
        let expected = sites_read_one_by_one("many.vcf", VcfOptions::default());
        assert_eq!(expected.len(), 475);
        for num_vars_per_block in [1, 7, 100, 1000] {
            let blocks = blocks_read(
                "many.vcf",
                VcfOptions::default(),
                Needs::CHROM_POS,
                Some(num_vars_per_block),
            );
            let given = sites_of(&blocks);
            assert_eq!(
                given.len(),
                expected.len(),
                "blocks of {num_vars_per_block} variants: the number of variants"
            );
            for (index, (given, expected)) in given.iter().zip(&expected).enumerate() {
                assert_eq!(
                    given, expected,
                    "blocks of {num_vars_per_block} variants: the variant {index}, counted from 0"
                );
            }
        }
    }

    #[test]
    fn a_collector_asked_for_the_genotypes_alone_gives_a_block_with_no_other_column() {
        // The VCF reader fills the chromosome and the position of every
        // variant it gives, and the block still has no column for them.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, Some(100));
        for block in &blocks {
            assert!(block.chrom.is_none());
            assert!(block.pos.is_none());
            assert!(block.id.is_none());
            assert!(block.alleles.is_none());
            assert!(block.qual.is_none());
        }
        assert_eq!(num_gts_of(&blocks), [10_000, 10_000, 10_000, 10_000, 7_500]);
    }

    #[test]
    fn the_default_number_of_variants_of_a_block_is_its_genotypes_between_the_two_bounds() {
        // 5 million genotypes divided by the individuals, which the
        // largest number of variants decides for 50 individuals and the
        // smallest one for 100000.
        assert_eq!(default_num_vars_per_block(50), 10_000);
        assert_eq!(default_num_vars_per_block(1000), 5_000);
        assert_eq!(default_num_vars_per_block(100_000), 100);
        // A source of no individual, which no VCF is: the division by one
        // individual gives the same, as pyNei's `max(num_samples, 1)`
        // does.
        assert_eq!(default_num_vars_per_block(0), 10_000);
        assert_eq!(default_num_vars_per_block(1), 10_000);
        assert_eq!(GENOTYPES_PER_BLOCK, 5_000_000);
        assert_eq!(MIN_NUM_VARS_PER_BLOCK, 100);
        assert_eq!(MAX_NUM_VARS_PER_BLOCK, 10_000);

        // The 50 individuals of `many.vcf` and no size asked for: its 475
        // variants are one block, since the default is 10000 of them.
        let blocks = blocks_read("many.vcf", VcfOptions::default(), Needs::GTS, None);
        assert_eq!(num_vars_of(&blocks), [475]);
    }

    /// The header of the VCFs written in these tests, with three
    /// individuals, and the first data line is its line 4.
    const HEADER: &str = "\
##fileformat=VCFv4.4
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3
";

    /// A VCF of the tests: the header above and one data line for each line
    /// given, whose columns are written with spaces and separated by tabs
    /// in the file.
    fn vcf_of(lines: &[&str]) -> String {
        let mut vcf = HEADER.to_string();
        for line in lines {
            vcf.push_str(&line.replace(' ', "\t"));
            vcf.push('\n');
        }
        vcf
    }

    /// A collector over a VCF written in a test.
    fn collector_over_text(
        vcf: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> BlockCollector<VcfReader<Cursor<Vec<u8>>>> {
        let reader = match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        match BlockCollector::new(reader, needs, num_vars_per_block) {
            Ok(collector) => collector,
            Err(error) => panic!("the collector was not built: {error}"),
        }
    }

    #[test]
    fn a_reader_with_no_variants_gives_no_block() {
        let mut collector =
            collector_over_text(HEADER, VcfOptions::default(), Needs::ALL, Some(100));
        assert!(collector.next_block().unwrap().is_none());
        assert!(collector.next_block().unwrap().is_none());
    }

    #[test]
    fn the_error_of_the_third_variant_comes_after_the_block_of_the_two_before_it() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 300 rs3 A T . PASS . GT 0/0 0/1 0/0/1/1",
        ]);
        let mut collector = collector_over_text(&vcf, VcfOptions::default(), Needs::GTS, Some(2));

        let first = collector.next_block().unwrap().expect("the first block");
        assert_eq!(first.num_vars, 2);
        assert_eq!(first.gts, [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]);

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::VcfGenotypePloidy {
            line,
            individual,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        // The third data line of a VCF of three header lines.
        assert_eq!(line, 6);
        assert_eq!(individual, "ind3");
        assert_eq!((found, expected), (4, 2));

        // The variant of the block that was being built is lost with the
        // error, and there is no block after it.
        assert!(collector.next_block().unwrap().is_none());
    }

    #[test]
    fn a_collector_of_blocks_of_no_variant_is_refused() {
        let reader = VcfReader::new(
            Cursor::new(HEADER.as_bytes().to_vec()),
            VcfOptions::default(),
        )
        .expect("the reader");
        let error = match BlockCollector::new(reader, Needs::GTS, Some(0)) {
            Ok(collector) => panic!("the collector was built: {collector:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");
        let message = error.to_string();
        assert!(message.contains('0'), "{message}");
    }

    /// The genotypes of a block are its variants times the individuals
    /// times the ploidy, which is a multiplication that a size asked for by
    /// a caller carries beyond what a `usize` holds. It is an error and not
    /// a panic, on a machine of 32 bit addresses as on one of 64.
    #[test]
    fn a_block_of_more_genotypes_than_a_usize_holds_is_refused_when_the_collector_is_built() {
        // Two individuals of the ploidy 2 are four genotypes in every
        // variant, so every size above a quarter of `usize::MAX` is one.
        let error = match BlockCollector::new(
            FakeReader::giving(0, Fills::Everything),
            Needs::GTS,
            Some(usize::MAX),
        ) {
            Ok(collector) => panic!("the collector was built: {collector:?}"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (num_vars_per_block, num_individuals, ploidy),
            (usize::MAX, 2, 2)
        );
    }

    /// The columns of a block are asked of the machine when the block is
    /// started, before a variant is read, and a size that a caller wrote
    /// reaches neither an abort nor a panic. A source of no individual
    /// holds no genotype, and its positions are 8 bytes in every variant,
    /// so the genotypes of its blocks are 0 and their memory is not.
    #[test]
    fn a_block_of_more_memory_than_the_machine_gives_is_refused_before_a_variant_is_read() {
        let mut collector = BlockCollector::new(
            FakeReader::of_no_individual(4),
            Needs::CHROM_POS,
            Some(usize::MAX),
        )
        .expect("the collector");

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (num_vars_per_block, num_individuals, ploidy),
            (usize::MAX, 0, 2)
        );

        assert_eq!(collector.reader().reads, 0);
        assert!(collector.next_block().expect("no block").is_none());
    }

    /// The alleles of a block are one buffer of text with the end of each
    /// allele, and the buffers are reserved when the block is started: a
    /// column that grows while it is filled is the allocations that one
    /// buffer is there to spare.
    #[test]
    fn the_buffers_of_an_alleles_column_hold_a_block_of_biallelic_variants() {
        let mut column = AllelesColumn::with_num_vars(1000).expect("the column");
        let texts = column.texts.capacity();
        let allele_ends = column.allele_ends.capacity();
        let var_ends = column.var_ends.capacity();

        let alleles = ["A".to_string(), "T".to_string()];
        for _ in 0..1000 {
            column.push(&alleles);
        }

        assert_eq!(column.texts.capacity(), texts);
        assert_eq!(column.allele_ends.capacity(), allele_ends);
        assert_eq!(column.var_ends.capacity(), var_ends);
        assert_eq!(column.num_vars(), 1000);
        assert_eq!(column.num_alleles(999), 2);
        assert_eq!(column.allele(999, 1), "T");
    }

    /// The quality of each variant of a block, with the NaN of a variant
    /// that has none as `None`, so that the assertion compares no float
    /// with another.
    fn quals_of(block: &Block) -> Vec<Option<f32>> {
        block
            .qual
            .as_ref()
            .expect("the qualities were asked for")
            .iter()
            .map(|qual| if qual.is_nan() { None } else { Some(*qual) })
            .collect()
    }

    /// The alleles of one variant of a block.
    fn alleles_of(block: &Block, var: usize) -> Vec<&str> {
        let alleles = block.alleles.as_ref().expect("the alleles were asked for");
        (0..alleles.num_alleles(var))
            .map(|allele| alleles.allele(var, allele))
            .collect()
    }

    #[test]
    fn every_column_of_cases_vcf_is_in_its_blocks() {
        // The four variants of the table of `cases.vcf` of "How it is
        // verified" of `docs/specs/io_vcf.md`, in blocks of three.
        let blocks = blocks_read("cases.vcf", every_variant(), Needs::ALL, Some(3));
        assert_eq!(num_vars_of(&blocks), [3, 1]);

        let first = &blocks[0];
        assert_eq!(first.num_individuals, 3);
        assert_eq!(first.ploidy, 2);
        // The chromosome of the four variants is `chr1`, the first name of
        // the file and so the number 0.
        assert_eq!(first.chrom.as_deref(), Some([0, 0, 0].as_slice()));
        assert_eq!(first.pos.as_deref(), Some([100, 200, 300].as_slice()));
        assert_eq!(
            first.id.as_deref(),
            Some(["rs1".to_string(), String::new(), String::new()].as_slice())
        );
        assert_eq!(quals_of(first), [Some(29.5), None, Some(67.0)]);
        assert_eq!(alleles_of(first, 0), ["A", "T"]);
        assert_eq!(alleles_of(first, 1), ["A", "T"]);
        assert_eq!(alleles_of(first, 2), ["A", "G", "T"]);
        // The genotypes of the three variants, `0/0 0/1 1/1`, `./. 0|1 .|0`
        // and `1/2 2|1 2/2`, one row for each of them.
        let gts: Vec<i8> = [
            [0, 0, 0, 1, 1, 1],
            [MISSING_ALLELE, MISSING_ALLELE, 0, 1, MISSING_ALLELE, 0],
            [1, 2, 2, 1, 2, 2],
        ]
        .concat();
        assert_eq!(first.gts, gts);

        let last = &blocks[1];
        assert_eq!(last.chrom.as_deref(), Some([0].as_slice()));
        assert_eq!(last.pos.as_deref(), Some([400].as_slice()));
        assert_eq!(last.id.as_deref(), Some([String::new()].as_slice()));
        assert_eq!(quals_of(last), [Some(47.0)]);
        // The variant with no alternative allele: the reference alone.
        assert_eq!(alleles_of(last, 0), ["T"]);
        assert_eq!(last.gts, [0, 0, 0, 0, 0, 0]);

        // A variant and an allele the column does not hold.
        let alleles = last.alleles.as_ref().expect("the alleles");
        assert_eq!(alleles.num_vars(), 1);
        assert_eq!(alleles.num_alleles(1), 0);
        assert_eq!(alleles.allele(0, 1), "");
        assert_eq!(alleles.allele(1, 0), "");
    }

    #[test]
    fn a_block_of_a_tetraploid_vcf_holds_four_alleles_for_each_individual() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0/1/1 0/1/1/1 ./././.",
            "chr1 200 rs2 A T . PASS . GT 0/0/0/0 1/1/1/1 0/0/0/1",
            "chr1 300 rs3 A T . PASS . GT 0/0/0/0 1/1/1/1 0/0/0/1",
        ]);
        let options = VcfOptions {
            ploidy: 4,
            only_passed: true,
        };
        let mut collector = collector_over_text(&vcf, options, Needs::GTS, Some(2));
        let blocks = blocks_of(&mut collector).expect("the blocks");

        assert_eq!(num_vars_of(&blocks), [2, 1]);
        assert_eq!(blocks[0].ploidy, 4);
        assert_eq!(blocks[0].num_individuals, 3);
        // 2 variants x 3 individuals x 4 alleles, and 1 x 3 x 4.
        assert_eq!(num_gts_of(&blocks), [24, 12]);
        assert_eq!(
            blocks[0].gts.get(..8),
            Some([0, 0, 1, 1, 0, 1, 1, 1].as_slice())
        );
        assert_eq!(
            blocks[0].gts.get(8..12),
            Some([MISSING_ALLELE; 4].as_slice())
        );
    }

    /// What a reader written for these tests puts in the variants it
    /// gives.
    #[derive(Clone, Copy)]
    enum Fills {
        /// The chromosome, the position and the genotypes of its two
        /// individuals.
        Everything,
        /// The chromosome and the position alone, which is a reader whose
        /// source has no genotype to give.
        NoGenotypes,
        /// One allele fewer than the two individuals of the ploidy 2 have,
        /// which is a reader with a defect.
        OneAlleleTooFew,
    }

    /// A reader of two individuals of the ploidy 2, written for these
    /// tests: a VCF cannot be made to do what a collector has to stand.
    struct FakeReader {
        individuals: Vec<String>,
        chroms: ChromTable,
        /// How many variants it still has to give.
        left: usize,
        fills: Fills,
        /// The variant it fails at, counted from 1.
        fails_at: Option<usize>,
        /// Which variant it is about to give, counted from 1.
        next_var: usize,
        /// How many times it was asked for a variant.
        reads: usize,
        /// What the collector asked this reader for.
        needs: Needs,
    }

    impl FakeReader {
        /// A reader of `num_vars` variants that fills each of them as
        /// `fills` says.
        fn giving(num_vars: usize, fills: Fills) -> FakeReader {
            FakeReader {
                individuals: vec!["ind1".to_string(), "ind2".to_string()],
                chroms: ChromTable::new(),
                left: num_vars,
                fills,
                fails_at: None,
                next_var: 1,
                reads: 0,
                needs: Needs::empty(),
            }
        }

        /// A reader whose source has no individual, and so no genotype: a
        /// block of it holds the columns of its variants and nothing else.
        fn of_no_individual(num_vars: usize) -> FakeReader {
            FakeReader {
                individuals: Vec::new(),
                ..FakeReader::giving(num_vars, Fills::NoGenotypes)
            }
        }

        /// A reader of `num_vars` variants whose variant number `fails_at`,
        /// counted from 1, is an error, and which gives the variants after
        /// it. The trait says that a reader ends at its error, so a
        /// collector that reads on after one reads what no reader of popnei
        /// gives: this one shows what it would do with it.
        fn failing_at(num_vars: usize, fails_at: usize) -> FakeReader {
            FakeReader {
                fails_at: Some(fails_at),
                ..FakeReader::giving(num_vars, Fills::Everything)
            }
        }
    }

    impl VariantReader for FakeReader {
        fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
            var.clear();
            self.reads = self.reads.saturating_add(1);
            let Some(left) = self.left.checked_sub(1) else {
                return Ok(false);
            };
            self.left = left;
            let this_var = self.next_var;
            self.next_var = this_var.saturating_add(1);
            if self.fails_at == Some(this_var) {
                return Err(Error::Io(std::io::Error::other(
                    "the reader of the tests failed",
                )));
            }
            var.chrom = self.chroms.intern("chr1");
            var.pos = 100;
            var.filled = Needs::CHROM_POS;
            match self.fills {
                Fills::Everything => {
                    var.gts.extend_from_slice(&[0, 1, 1, 1]);
                    var.filled |= Needs::GTS;
                }
                Fills::NoGenotypes => {}
                Fills::OneAlleleTooFew => {
                    var.gts.extend_from_slice(&[0, 1, 1]);
                    var.filled |= Needs::GTS;
                }
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

    /// A column nobody wants is never parsed, so what the collector asks
    /// its reader for is what its blocks will hold.
    #[test]
    fn a_collector_asks_its_reader_for_its_columns_and_the_genotypes() {
        let no_variant = || FakeReader::giving(0, Fills::Everything);
        let collector =
            BlockCollector::new(no_variant(), Needs::QUAL, Some(2)).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::GTS | Needs::QUAL);

        let collector =
            BlockCollector::new(no_variant(), Needs::empty(), None).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::GTS);

        let collector =
            BlockCollector::new(no_variant(), Needs::ALL, Some(2)).expect("the collector");
        assert_eq!(collector.reader().needs, Needs::ALL);
    }

    /// The names a Python and a TypeScript user writes in `fields` are the
    /// names of the columns of a block, and which field of the core each
    /// one asks for is knowledge of the domain that the core keeps: the
    /// two binding crates call this with what their user gave.
    #[test]
    fn the_names_of_the_columns_of_a_block_ask_for_the_fields_of_the_core() {
        assert_eq!(FIELD_NAMES, ["chrom", "pos", "id", "alleles", "qual"]);
        assert_eq!(FIELD_NAMES, FIELDS_OF_THE_NAMES.map(|(name, _)| name));

        // Every block holds the genotypes, so they are part of what any
        // set of names asks for, and no name asks for them.
        let no_name: [&str; 0] = [];
        assert_eq!(needs_of_the_fields(no_name).expect("no name"), Needs::GTS);
        assert_eq!(
            needs_of_the_fields(["chrom"]).expect("the chromosomes"),
            Needs::GTS | Needs::CHROM_POS
        );
        // The chromosome and the position are one field of the core, so
        // either name fills both.
        assert_eq!(
            needs_of_the_fields(["pos"]).expect("the positions"),
            Needs::GTS | Needs::CHROM_POS
        );
        assert_eq!(
            needs_of_the_fields(["qual", "id", "alleles", "chrom"]).expect("the four"),
            Needs::ALL
        );
        assert_eq!(
            needs_of_the_fields(["id", "id"]).expect("the ids twice"),
            Needs::GTS | Needs::ID
        );

        let error = match needs_of_the_fields(["chrom", "c"]) {
            Ok(needs) => panic!("the names gave {needs}"),
            Err(error) => error,
        };
        let Error::NotAFieldOfABlock { name } = &error else {
            panic!("the error is {error}");
        };
        assert_eq!(name, "c");
        let message = error.to_string();
        assert!(
            message.contains("`c` is not a field of a block"),
            "{message}"
        );
        assert!(
            message.contains("`chrom`, `pos`, `id`, `alleles`, `qual`"),
            "{message}"
        );

        // The genotypes are not a name a user writes.
        let error = needs_of_the_fields(["gts"]);
        assert!(
            matches!(error, Err(Error::NotAFieldOfABlock { .. })),
            "`gts` gave {error:?}"
        );
    }

    /// The genotypes of a block are read as variants x individuals x
    /// ploidy, so a variant of another number of alleles would be read
    /// wrong: numpy refuses the reshape in Python, and in TypeScript the
    /// genotypes cross flat and nothing says that they moved.
    #[test]
    fn a_variant_of_a_number_of_alleles_other_than_the_individuals_is_an_error() {
        let mut collector = BlockCollector::new(
            FakeReader::giving(3, Fills::OneAlleleTooFew),
            Needs::empty(),
            Some(2),
        )
        .expect("the collector");

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::VariantOfAnotherSize {
            found,
            expected,
            num_individuals,
            ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((found, expected), (3, 4));
        assert_eq!((num_individuals, ploidy), (2, 2));

        assert!(collector.next_block().expect("no block").is_none());
    }

    /// The trait says that a reader ends at its error, and the collector
    /// does not lean on it: a reader that goes on giving variants after one
    /// gives the collector's caller no block after it.
    #[test]
    fn the_collector_gives_no_block_after_the_error_of_its_reader() {
        let mut collector =
            BlockCollector::new(FakeReader::failing_at(6, 3), Needs::CHROM_POS, Some(2))
                .expect("the collector");

        let first = collector
            .next_block()
            .expect("the first block")
            .expect("a block");
        assert_eq!(first.num_vars, 2);

        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::Io(_)), "the error is {error}");

        // The reader has three variants left and gives them, and the
        // collector gives no block.
        assert!(collector.next_block().expect("no block").is_none());
        assert!(collector.next_block().expect("no block").is_none());
    }

    /// One row of the table of `cases.vcf` of "How it is verified" of
    /// `docs/specs/io_vcf.md`, which the blocks of these tests are built
    /// from: three individuals of the ploidy 2, and the chromosome `chr1`,
    /// the first name of that file and so the number 0.
    struct Row {
        pos: u64,
        id: &'static str,
        alleles: &'static [&'static str],
        /// `None` for the variant whose QUAL is `.`, which is a NaN in the
        /// column of a block.
        qual: Option<f32>,
        gts: [i8; 6],
    }

    const MISSING: i8 = MISSING_ALLELE;

    /// The four variants of `cases.vcf`, in the order of the file.
    const CASES: [Row; 4] = [
        Row {
            pos: 100,
            id: "rs1",
            alleles: &["A", "T"],
            qual: Some(29.5),
            gts: [0, 0, 0, 1, 1, 1],
        },
        Row {
            pos: 200,
            id: "",
            alleles: &["A", "T"],
            qual: None,
            gts: [MISSING, MISSING, 0, 1, MISSING, 0],
        },
        Row {
            pos: 300,
            id: "",
            alleles: &["A", "G", "T"],
            qual: Some(67.0),
            gts: [1, 2, 2, 1, 2, 2],
        },
        Row {
            pos: 400,
            id: "",
            alleles: &["T"],
            qual: Some(47.0),
            gts: [0, 0, 0, 0, 0, 0],
        },
    ];

    /// A block built by hand from the rows of that table, with every
    /// column: this is what a reader of `cases.vcf` gives, and building it
    /// here is what lets the tests of the block and of `reblock` stand on
    /// the literals of the spec without a file.
    fn cases_block(rows: &[usize]) -> Block {
        let mut gts = Vec::new();
        let mut chrom = Vec::new();
        let mut pos = Vec::new();
        let mut id = Vec::new();
        let mut qual = Vec::new();
        let mut alleles = AllelesColumn::with_num_vars(rows.len()).expect("the alleles");
        for index in rows {
            let row = &CASES[*index];
            gts.extend_from_slice(&row.gts);
            chrom.push(0);
            pos.push(row.pos);
            id.push(row.id.to_string());
            qual.push(row.qual.unwrap_or(f32::NAN));
            let texts: Vec<String> = row.alleles.iter().map(|text| text.to_string()).collect();
            alleles.push(&texts);
        }
        Block {
            num_vars: rows.len(),
            num_individuals: 3,
            ploidy: 2,
            gts,
            chrom: Some(chrom),
            pos: Some(pos),
            id: Some(id),
            alleles: Some(alleles),
            qual: Some(qual),
        }
    }

    /// The alleles of a view, as texts.
    fn alleles_of_view<'a>(view: &VariantRef<'a>) -> Vec<&'a str> {
        let num_alleles = view.num_alleles().expect("the alleles were asked for");
        (0..num_alleles)
            .map(|allele| view.allele(allele).unwrap_or(""))
            .collect()
    }

    /// That a view holds the row `row` of the table of `cases.vcf`, every
    /// field of it. The quality is compared as an `Option`, with the NaN
    /// of a variant that has none as `None`.
    fn assert_view_is_the_row(view: &VariantRef<'_>, row: usize) {
        let expected = &CASES[row];
        assert_eq!(view.gts(), expected.gts, "the genotypes of the row {row}");
        assert_eq!(view.chrom(), Some(0), "the chromosome of the row {row}");
        assert_eq!(
            view.pos(),
            Some(expected.pos),
            "the position of the row {row}"
        );
        assert_eq!(view.id(), Some(expected.id), "the id of the row {row}");
        let qual = view.qual().filter(|qual| !qual.is_nan());
        assert_eq!(qual, expected.qual, "the quality of the row {row}");
        assert_eq!(
            alleles_of_view(view),
            expected.alleles,
            "the alleles of the row {row}"
        );
    }

    /// The views of a block are its variants in order, each with every
    /// field of its row.
    #[test]
    fn the_views_of_a_block_give_the_fields_of_each_of_its_variants() {
        let block = cases_block(&[0, 1, 2, 3]);
        assert_eq!(block.fields(), Needs::ALL);
        block.check().expect("the block is of its size");

        let views: Vec<VariantRef<'_>> = block.variants().collect();
        assert_eq!(views.len(), 4);
        for (row, view) in views.iter().enumerate() {
            assert_view_is_the_row(view, row);
        }

        // The view of one variant, and no view for a fifth variant of a
        // block of four.
        let third = block.variant(2).expect("the third variant");
        assert_view_is_the_row(&third, 2);
        assert!(block.variant(4).is_none());

        // An allele the variant does not have: the fourth row has the
        // reference allele alone.
        let last = block.variant(3).expect("the fourth variant");
        assert_eq!(last.num_alleles(), Some(1));
        assert_eq!(last.allele(0), Some("T"));
        assert_eq!(last.allele(1), None);
    }

    /// A column that is not in the block is `None` in every view and in
    /// the fields of the block, which is what a consumer that depends on a
    /// field reads.
    #[test]
    fn the_views_of_a_block_give_none_for_a_column_it_does_not_hold() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        block.id = None;
        block.qual = None;
        block.alleles = None;
        assert_eq!(block.fields(), Needs::GTS | Needs::CHROM_POS);
        assert_eq!(
            Needs::ALL.difference(block.fields()),
            Needs::ID | Needs::ALLELES | Needs::QUAL
        );

        let first = block.variant(0).expect("the first variant");
        assert_eq!(first.gts(), CASES[0].gts);
        assert_eq!(first.pos(), Some(100));
        assert_eq!(first.id(), None);
        assert_eq!(first.qual(), None);
        assert_eq!(first.num_alleles(), None);
        assert_eq!(first.allele(0), None);

        // A block built without the genotypes: its views have none, and
        // its variants are still counted.
        block.gts.clear();
        assert_eq!(block.fields(), Needs::CHROM_POS);
        let first = block.variant(0).expect("the first variant");
        assert!(first.gts().is_empty());
        assert_eq!(first.pos(), Some(100));
        assert_eq!(block.variants().count(), 4);
    }

    /// The genotypes of a block are read as variants x individuals x
    /// ploidy, so a block whose array is not of that size would be read
    /// with every genotype after the fault at the place of another.
    #[test]
    fn a_block_whose_genotypes_lost_an_allele_fails_check() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        block.gts.pop();
        let error = match block.check() {
            Ok(()) => panic!("the block passed the check"),
            Err(error) => error,
        };
        let Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        // 4 variants x 3 individuals x 2 alleles.
        assert_eq!((array, found, expected), ("gts", 23, 24));
        // A user who gets this reads which array and both sizes.
        let message = error.to_string();
        assert!(message.contains("`gts`"), "{message}");
        assert!(message.contains("23"), "{message}");
        assert!(message.contains("24"), "{message}");

        // A column of another size than the variants of the block.
        let mut block = cases_block(&[0, 1, 2, 3]);
        block.pos.as_mut().expect("the positions").pop();
        let error = match block.check() {
            Ok(()) => panic!("the block passed the check"),
            Err(error) => error,
        };
        let Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((array, found, expected), ("pos", 3, 4));
    }

    /// A filter keeps some variants of a block and the block is compacted
    /// in place: every column keeps the rows that stayed, in their order.
    #[test]
    fn keeping_three_of_the_four_variants_leaves_a_block_of_those_three() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        let capacity = block.gts.capacity();
        block
            .retain_vars(&[true, false, true, true])
            .expect("the variants to keep");

        assert_eq!(block.num_vars, 3);
        block.check().expect("the block is of its size");
        assert_eq!(block.gts.capacity(), capacity);
        let views: Vec<VariantRef<'_>> = block.variants().collect();
        assert_eq!(views.len(), 3);
        // The rows 1, 3 and 4 of the table, counted from 1.
        for (view, row) in views.iter().zip([0, 2, 3]) {
            assert_view_is_the_row(view, row);
        }
    }

    /// A filter that drops every variant of a block: the block is left
    /// empty and its reader gives it to nobody.
    #[test]
    fn keeping_no_variant_leaves_a_block_of_no_variant_and_empty_columns() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        block
            .retain_vars(&[false, false, false, false])
            .expect("the variants to keep");

        assert_eq!(block.num_vars, 0);
        block.check().expect("the block is of its size");
        assert!(block.gts.is_empty());
        assert_eq!(block.chrom.as_deref(), Some([].as_slice()));
        assert_eq!(block.pos.as_deref(), Some([].as_slice()));
        assert_eq!(block.id.as_deref(), Some([].as_slice()));
        assert_eq!(block.qual.as_ref().map(Vec::len), Some(0));
        assert_eq!(block.alleles.as_ref().map(AllelesColumn::num_vars), Some(0));
        assert_eq!(block.variants().count(), 0);
    }

    /// A `keep` of another size is a defect of the filter that wrote it,
    /// and the block is not touched: the variants it holds are still the
    /// ones its reader gave.
    #[test]
    fn a_keep_of_three_values_for_four_variants_is_an_error_and_leaves_the_block_as_it_was() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        let error = match block.retain_vars(&[true, false, true]) {
            Ok(()) => panic!("the block kept 3 values for 4 variants"),
            Err(error) => error,
        };
        let Error::KeepOfAnotherSize { found, num_vars } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((found, num_vars), (3, 4));

        assert_eq!(block.num_vars, 4);
        let views: Vec<VariantRef<'_>> = block.variants().collect();
        assert_eq!(views.len(), 4);
        for (row, view) in views.iter().enumerate() {
            assert_view_is_the_row(view, row);
        }
    }

    /// A reader of blocks written for these tests: it gives the blocks it
    /// was built with, keeps the address of the genotypes of each one, so
    /// that a test sees whether a block was copied, and counts the calls,
    /// so that a test sees that a source is not called again after its
    /// error.
    struct GivenBlocks {
        individuals: Vec<String>,
        chroms: ChromTable,
        ploidy: usize,
        /// The blocks still to give, the last one first.
        left: Vec<Block>,
        /// Where the genotypes of every block it gave are in memory.
        addresses: Vec<usize>,
        /// The call at which it gives an error instead of a block, counted
        /// from 1. It has blocks after it, which the trait says no reader
        /// of popnei gives; this one gives them to show what `reblock`
        /// does with a source that goes on.
        fails_at: Option<usize>,
        /// How many times it was asked for a block.
        calls: usize,
        /// What it was last asked to fill.
        needs: Needs,
    }

    impl GivenBlocks {
        /// A reader of three individuals of the ploidy 2, which is what
        /// `cases.vcf` has, that gives `blocks` in their order.
        fn of(blocks: Vec<Block>) -> GivenBlocks {
            let mut chroms = ChromTable::new();
            // The chromosome of the four rows of `cases.vcf`, whose number
            // is 0 in the blocks built from them.
            chroms.intern("chr1");
            let mut left = blocks;
            left.reverse();
            GivenBlocks {
                individuals: vec!["ind1".to_string(), "ind2".to_string(), "ind3".to_string()],
                chroms,
                ploidy: 2,
                left,
                addresses: Vec::new(),
                fails_at: None,
                calls: 0,
                needs: Needs::ALL,
            }
        }

        /// The same reader, whose call number `call` is an error.
        fn failing_at(blocks: Vec<Block>, call: usize) -> GivenBlocks {
            GivenBlocks {
                fails_at: Some(call),
                ..GivenBlocks::of(blocks)
            }
        }

        /// One block for each row of the table of `cases.vcf`.
        fn of_one_variant_each() -> GivenBlocks {
            GivenBlocks::of(vec![
                cases_block(&[0]),
                cases_block(&[1]),
                cases_block(&[2]),
                cases_block(&[3]),
            ])
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> Result<Option<Block>> {
            self.calls = self.calls.saturating_add(1);
            if self.fails_at == Some(self.calls) {
                return Err(Error::Io(std::io::Error::other(
                    "the reader of the tests failed",
                )));
            }
            let Some(block) = self.left.pop() else {
                return Ok(None);
            };
            self.addresses.push(block.gts.as_ptr().addr());
            Ok(Some(block))
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
    }

    /// Every block a reader gives, until it has no more or it fails.
    fn blocks_given(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// `reblock` joins the blocks that are too short: four blocks of one
    /// variant, in blocks of three, are one block of three variants and
    /// one of one, with every column of the four rows.
    #[test]
    fn reblock_of_three_joins_four_blocks_of_one_variant_into_one_of_three_and_one_of_one() {
        let mut reblock =
            Reblock::new(GivenBlocks::of_one_variant_each(), Some(3)).expect("the reblock");
        let blocks = blocks_given(&mut reblock).expect("the blocks");

        assert_eq!(num_vars_of(&blocks), [3, 1]);
        let first: Vec<VariantRef<'_>> = blocks[0].variants().collect();
        assert_eq!(first.len(), 3);
        for (row, view) in first.iter().enumerate() {
            assert_view_is_the_row(view, row);
        }
        blocks[0].check().expect("the joined block is of its size");
        let last: Vec<VariantRef<'_>> = blocks[1].variants().collect();
        assert_eq!(last.len(), 1);
        assert_view_is_the_row(&last[0], 3);
        blocks[1].check().expect("the last block is of its size");

        // The source gave its four blocks and was asked once more, which
        // is the call that told `reblock` there were no more.
        assert_eq!(reblock.reader.calls, 5);
        assert!(reblock.next_block().expect("no block").is_none());
    }

    /// A block that already has the size asked for, with nothing waiting
    /// from the block before it, goes through as it is: its genotypes are
    /// at the address the source allocated, and no row was copied.
    #[test]
    fn reblock_of_one_gives_the_blocks_of_one_variant_of_its_source_with_no_copy() {
        let mut reblock =
            Reblock::new(GivenBlocks::of_one_variant_each(), Some(1)).expect("the reblock");
        let blocks = blocks_given(&mut reblock).expect("the blocks");

        assert_eq!(num_vars_of(&blocks), [1, 1, 1, 1]);
        let given: Vec<usize> = blocks
            .iter()
            .map(|block| block.gts.as_ptr().addr())
            .collect();
        assert_eq!(given, reblock.reader.addresses);
        for (row, block) in blocks.iter().enumerate() {
            let view = block.variant(0).expect("the variant of the block");
            assert_view_is_the_row(&view, row);
        }
    }

    /// The trait says that a reader ends at its error and that a reader
    /// over another reader does not call its source again: what a source
    /// gives after an error is not the file.
    #[test]
    fn reblock_calls_a_source_that_failed_once_and_no_more() {
        let blocks = vec![cases_block(&[0]), cases_block(&[1]), cases_block(&[2])];
        let mut reblock =
            Reblock::new(GivenBlocks::failing_at(blocks, 2), Some(3)).expect("the reblock");

        let error = match reblock.next_block() {
            Ok(block) => panic!("the reblock gave {block:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::Io(_)), "the error is {error}");
        // The first block was waiting to be joined and is lost with the
        // error, as the variants that do not fill a block are.
        assert_eq!(reblock.reader.calls, 2);
        assert!(reblock.next_block().expect("no block").is_none());
        assert!(reblock.next_block().expect("no block").is_none());
        assert_eq!(reblock.reader.calls, 2);
    }

    /// The rows of two blocks are joined into one array of variants x
    /// individuals x ploidy, so blocks of another number of individuals or
    /// another ploidy cannot be joined.
    #[test]
    fn reblock_refuses_blocks_that_do_not_fit_together() {
        let mut two_individuals = cases_block(&[1]);
        two_individuals.num_individuals = 2;
        two_individuals.gts.truncate(4);
        let blocks = vec![cases_block(&[0]), two_individuals];
        let mut reblock = Reblock::new(GivenBlocks::of(blocks), Some(3)).expect("the reblock");

        let error = match reblock.next_block() {
            Ok(block) => panic!("the reblock gave {block:?}"),
            Err(error) => error,
        };
        let Error::BlocksDoNotFitTogether {
            num_individuals,
            ploidy,
            found_num_individuals,
            found_ploidy,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((num_individuals, ploidy), (3, 2));
        assert_eq!((found_num_individuals, found_ploidy), (2, 2));
        let message = error.to_string();
        assert!(message.contains("2 individuals"), "{message}");
        assert!(message.contains("3 individuals"), "{message}");
        assert!(reblock.next_block().expect("no block").is_none());
    }

    /// `reblock` checks every block it takes, because it joins their rows
    /// into one array and a block whose arrays are not of its size would
    /// move every row after it.
    #[test]
    fn reblock_refuses_a_block_whose_arrays_are_not_of_its_size() {
        let mut short = cases_block(&[0]);
        short.gts.pop();
        let mut reblock = Reblock::new(GivenBlocks::of(vec![short]), Some(1)).expect("the reblock");

        let error = match reblock.next_block() {
            Ok(block) => panic!("the reblock gave {block:?}"),
            Err(error) => error,
        };
        let Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((array, found, expected), ("gts", 5, 6));
        assert!(reblock.next_block().expect("no block").is_none());
    }

    /// A change of `Needs` in the middle of a pass changes the columns of
    /// the blocks that come after it, and blocks of different columns are
    /// not joined: what was waiting is given as a shorter block.
    #[test]
    fn reblock_gives_a_shorter_block_when_the_columns_of_the_source_change() {
        let mut with_the_genotypes_alone = cases_block(&[1]);
        with_the_genotypes_alone.chrom = None;
        with_the_genotypes_alone.pos = None;
        with_the_genotypes_alone.id = None;
        with_the_genotypes_alone.alleles = None;
        with_the_genotypes_alone.qual = None;
        let blocks = vec![cases_block(&[0]), with_the_genotypes_alone];
        let mut reblock = Reblock::new(GivenBlocks::of(blocks), Some(3)).expect("the reblock");

        // What the caller asks for reaches the source.
        reblock.set_needs(Needs::GTS);
        assert_eq!(reblock.reader.needs, Needs::GTS);

        let given = blocks_given(&mut reblock).expect("the blocks");
        assert_eq!(num_vars_of(&given), [1, 1]);
        assert_eq!(given[0].fields(), Needs::ALL);
        assert_view_is_the_row(&given[0].variant(0).expect("the variant"), 0);
        assert_eq!(given[1].fields(), Needs::GTS);
        assert_eq!(
            given[1].variant(0).expect("the variant").gts(),
            CASES[1].gts
        );
    }

    /// Both binding crates hold a `Box<dyn BlockReader>` and put a
    /// `Reblock` over it, so one is read here.
    #[test]
    fn a_boxed_reader_of_blocks_is_read_through_reblock() {
        let reader: Box<dyn BlockReader> = Box::new(GivenBlocks::of_one_variant_each());
        let mut reblock = Reblock::new(reader, Some(2)).expect("the reblock");

        assert_eq!(reblock.individuals(), ["ind1", "ind2", "ind3"]);
        assert_eq!(reblock.ploidy(), 2);
        assert_eq!(reblock.chroms().name(0), Some("chr1"));

        let blocks = blocks_given(&mut reblock).expect("the blocks");
        assert_eq!(num_vars_of(&blocks), [2, 2]);
        assert_view_is_the_row(&blocks[0].variant(0).expect("the variant"), 0);
        assert_view_is_the_row(&blocks[0].variant(1).expect("the variant"), 1);
        assert_view_is_the_row(&blocks[1].variant(0).expect("the variant"), 2);
        assert_view_is_the_row(&blocks[1].variant(1).expect("the variant"), 3);
    }

    /// Every reader keeps the rule that a block it gives holds one variant
    /// at least, and `reblock` does not trust its source to keep it: a
    /// source that gives a block of no variants has a defect, and asking
    /// it again would spin for ever over a source that always gives one.
    #[test]
    fn reblock_refuses_a_block_of_no_variants_and_asks_its_source_no_more() {
        // The block of no variants comes first and a block of one variant
        // after it, which a `reblock` that went on would give.
        let blocks = vec![cases_block(&[]), cases_block(&[0])];
        let mut reblock = Reblock::new(GivenBlocks::of(blocks), Some(1)).expect("the reblock");

        let error = match reblock.next_block() {
            Ok(block) => panic!("the reblock gave {block:?}"),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::ReaderGaveABlockOfNoVariants),
            "the error is {error}"
        );
        let message = error.to_string();
        assert!(message.contains("no variants"), "{message}");

        assert_eq!(reblock.reader.calls, 1);
        assert!(reblock.next_block().expect("no block").is_none());
        assert_eq!(reblock.reader.calls, 1);
    }

    /// A block holds one variant at least, so a reader that was asked for
    /// blocks of 0 variants is refused where it is built.
    #[test]
    fn a_reblock_of_blocks_of_no_variant_is_refused() {
        let error = match Reblock::new(GivenBlocks::of(Vec::new()), Some(0)) {
            Ok(reblock) => panic!("the reblock was built: {reblock:?}"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");
        let message = error.to_string();
        assert!(message.contains('0'), "{message}");
    }

    /// The blocks of one of the reference VCFs, read through the collector
    /// that exists, which is a reader of blocks.
    fn collected_over(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> CollectedBlocks<VcfReader<BufReader<File>>> {
        let reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        match CollectedBlocks::new(reader, needs, num_vars_per_block) {
            Ok(blocks) => blocks,
            Err(error) => panic!("{name}: the reader of blocks was not built: {error}"),
        }
    }

    /// `reblock` over a reader of a real file: the 500 variants of
    /// `many.vcf`, given in blocks of 100, come out in the four sizes of
    /// "How it is verified" of `docs/specs/block.md`, and what they hold,
    /// joined, is the same at the four.
    #[test]
    fn reblock_gives_the_variants_of_many_vcf_in_the_size_that_was_asked_for() {
        let mut expected_of_seven = vec![7; 71];
        expected_of_seven.push(3);
        let sizes: [(Option<usize>, Vec<usize>); 5] = [
            (Some(7), expected_of_seven),
            (Some(1), vec![1; 500]),
            (Some(1000), vec![500]),
            (Some(100), vec![100; 5]),
            // No size asked for: 10000 variants for 50 individuals, so the
            // 500 of the file are one block.
            (None, vec![500]),
        ];
        let mut first_sites: Option<Vec<Site>> = None;
        for (num_vars_per_block, expected) in sizes {
            let source = collected_over("many.vcf", every_variant(), Needs::CHROM_POS, Some(100));
            let mut reblock = Reblock::new(source, num_vars_per_block).expect("the reblock");
            let blocks = match blocks_given(&mut reblock) {
                Ok(blocks) => blocks,
                Err(error) => panic!("blocks of {num_vars_per_block:?}: {error}"),
            };
            assert_eq!(
                num_vars_of(&blocks),
                expected,
                "the variants of the blocks of {num_vars_per_block:?}"
            );
            // The genotypes, the positions and the chromosome numbers of
            // the blocks, joined, one variant after another.
            let sites = sites_of(&blocks);
            assert_eq!(sites.len(), 500);
            match first_sites.as_ref() {
                None => first_sites = Some(sites),
                Some(first) => assert_eq!(
                    &sites, first,
                    "blocks of {num_vars_per_block:?}: the variants are not the ones of the blocks of 7"
                ),
            }
        }
    }

    /// The two binding crates hold their collector over a boxed reader,
    /// `BlockCollector<Box<dyn VariantReader>>`, so one is built here.
    #[test]
    fn a_field_that_was_asked_for_and_that_the_reader_does_not_fill_is_an_error() {
        let reader: Box<dyn VariantReader> = Box::new(FakeReader::giving(2, Fills::NoGenotypes));
        let mut collector =
            BlockCollector::new(reader, Needs::CHROM_POS, Some(2)).expect("the collector");
        assert_eq!(collector.reader().individuals(), ["ind1", "ind2"]);
        assert_eq!(collector.reader().ploidy(), 2);
        let error = match collector.next_block() {
            Ok(block) => panic!("the collector gave {block:?}"),
            Err(error) => error,
        };
        let Error::FieldsNotFilled { fields } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(fields, Needs::GTS);
    }
}
