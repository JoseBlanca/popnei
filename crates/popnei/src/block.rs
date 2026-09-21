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
//! `docs/specs/block.md` has the design and section 2 of
//! `docs/architecture.md` the reasons for the arrays.

use std::collections::TryReserveError;
use std::fmt;

use crate::error::{Error, Result};
use crate::filters::FilteringStats;
use crate::variant::{ChromTable, Needs, VariantRef};

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

/// Which of the three sizes of a block a reader is working with.
///
/// What a caller does about a block the machine cannot give the memory for
/// depends on it, so the error carries it: a caller who asked for a size
/// asks for fewer variants, one who asked for none learns that the size
/// popnei chose for these individuals does not fit and passes one that
/// does, and one reading a file whose batches fix the size writes that file
/// again with smaller ones.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockSize {
    /// The `num_vars_per_block` that the caller wrote.
    AskedFor,
    /// [`default_num_vars_per_block`] for the individuals of the source,
    /// which is what a caller who asked for no size gets.
    ChosenByPopnei,
    /// The size a file fixed: the batch of a vars file, which its reader
    /// builds whole whatever size the caller asked its blocks to be.
    FixedByAFile,
}

impl BlockSize {
    /// What a caller does about a block of this size that does not fit,
    /// which is the end of the message of [`Error::BlockTooLarge`].
    pub(crate) fn way_out(self) -> &'static str {
        match self {
            BlockSize::AskedFor => "ask for fewer variants in a block",
            BlockSize::ChosenByPopnei => {
                "popnei chose that size for these individuals and this ploidy; \
                 ask for the blocks with a `num_vars_per_block` that fits"
            }
            BlockSize::FixedByAFile => {
                "the batches of the file fix that size, whatever size its blocks are asked \
                 for; write the file again with a smaller `num_vars_per_block`"
            }
        }
    }
}

/// How many variants a block of a reader over a source of `num_individuals`
/// individuals of `ploidy` holds: the size the caller asked for, or
/// [`default_num_vars_per_block`] when they asked for none, and which of
/// the two it is.
///
/// Every reader that takes a size makes this check, and nothing of a source
/// is read before it: a size that a caller wrote reaches neither an abort
/// nor a panic.
///
/// # Errors
///
/// When the size is 0, and when the genotypes of one block are more than a
/// `usize` counts.
pub(crate) fn size_of_the_blocks(
    num_vars_per_block: Option<usize>,
    num_individuals: usize,
    ploidy: usize,
) -> Result<(usize, BlockSize)> {
    let (num_vars_per_block, size) = match num_vars_per_block {
        Some(asked_for) => (asked_for, BlockSize::AskedFor),
        None => (
            default_num_vars_per_block(num_individuals),
            BlockSize::ChosenByPopnei,
        ),
    };
    check_the_size_of_a_block(num_vars_per_block, num_individuals, ploidy, size)?;
    Ok((num_vars_per_block, size))
}

/// The genotypes of a full block, `num_vars_per_block` times
/// `num_individuals` times `ploidy`, when the machine counts that many.
///
/// # Errors
///
/// When `num_vars_per_block` is 0, which is no block at all, and when the
/// multiplication carries beyond what a `usize` holds, which in wasm, where
/// a `usize` is 32 bits, is 10000 variants of 250000 individuals of the
/// ploidy 2.
pub(crate) fn check_the_size_of_a_block(
    num_vars_per_block: usize,
    num_individuals: usize,
    ploidy: usize,
    size: BlockSize,
) -> Result<usize> {
    if num_vars_per_block == 0 {
        return Err(Error::BlockOfNoVariants);
    }
    num_individuals
        .checked_mul(ploidy)
        .and_then(|alleles_per_var| alleles_per_var.checked_mul(num_vars_per_block))
        .ok_or(Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
            size,
        })
}

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

/// The fields a reader has to be asked for to fill the columns that
/// `names` name, the genotypes among them, which every block a user gets
/// holds.
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
    /// Every reader of the crate that fills the alleles of a block starts
    /// its column here, with the variants a full block holds, and then
    /// pushes the alleles of one variant after another with
    /// [`AllelesColumn::push`]. The memory is asked for with `try_reserve`,
    /// which gives it back as an error: `Vec::with_capacity` ends the
    /// process when the machine has not the memory, and a size that a
    /// caller of popnei wrote reaches it.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the buffers, which the
    /// reader turns into the error of the crate that names the size that
    /// was asked for.
    pub(crate) fn with_num_vars(
        num_vars: usize,
    ) -> std::result::Result<AllelesColumn, TryReserveError> {
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

    /// The alleles of one more variant, at the end of the column, the
    /// reference allele first and then the alternative ones, as the source
    /// gives their texts. It is how a reader fills a column, one variant
    /// after another and in the order of the block.
    ///
    /// The buffers grow when the alleles of a block are more, or longer,
    /// than [`AllelesColumn::with_num_vars`] reserved for, and growing is
    /// what `Vec` does with a memory it is not given: an allocation that
    /// fails here ends the process. A reader that reserves for the block it
    /// builds does not reach it.
    pub(crate) fn push(&mut self, alleles: &[String]) {
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

    /// Where the alleles of the variants before `var` end in `allele_ends`:
    /// the first allele of the variant `var`.
    fn allele_before(&self, var: usize) -> usize {
        var.checked_sub(1)
            .and_then(|before| self.var_ends.get(before).copied())
            .unwrap_or(0)
    }

    /// The alleles of `count` variants from `from` on, copied into a column
    /// of their own, which is what cutting a block does with them. This
    /// column is left as it is.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the new column.
    fn try_rows(
        &self,
        from: usize,
        count: usize,
    ) -> std::result::Result<AllelesColumn, TryReserveError> {
        let from = from.min(self.var_ends.len());
        let end = from.saturating_add(count).min(self.var_ends.len());
        let allele_from = self.allele_before(from);
        let allele_end = self.allele_before(end).max(allele_from);
        let byte_from = self.byte_before(from);
        let byte_end = allele_end
            .checked_sub(1)
            .and_then(|last| self.allele_ends.get(last).copied())
            .unwrap_or(byte_from);
        let mut texts = Vec::new();
        let bytes = self.texts.get(byte_from..byte_end).unwrap_or(&[]);
        texts.try_reserve_exact(bytes.len())?;
        texts.extend_from_slice(bytes);
        let ends = self.allele_ends.get(allele_from..allele_end).unwrap_or(&[]);
        let mut allele_ends = Vec::new();
        allele_ends.try_reserve_exact(ends.len())?;
        allele_ends.extend(ends.iter().map(|end| end.saturating_sub(byte_from)));
        let ends = self.var_ends.get(from..end).unwrap_or(&[]);
        let mut var_ends = Vec::new();
        var_ends.try_reserve_exact(ends.len())?;
        var_ends.extend(ends.iter().map(|end| end.saturating_sub(allele_from)));
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
/// A column other than the genotypes is there only when the reader was
/// asked for it, and `None` when it was not or when its source has no such
/// field. The number of a chromosome is a number of the table of the reader
/// the block came from, and that table grows while the source is read, so
/// the name of a number is looked up after the block was given.
///
/// The fields are public, because every reader builds blocks, so a reader
/// with a defect can build one whose arrays are not of its size.
/// [`Block::check`] is what says whether a block keeps them.
#[derive(Debug)]
pub struct Block {
    /// How many variants the block holds. A reader that was given a size
    /// gives that many, the last block of a source aside; a filter and the
    /// reader of a vars file give the blocks the size their work leaves.
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
    /// [`ChromTable`] of the reader the block
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
        // Every field of the block is named here, and in the five other
        // places that read them all, so that a column added later does not
        // fall out of one of them without the compiler saying so.
        let Block {
            num_vars,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = self;
        let mut fields = Needs::empty();
        if !gts.is_empty() || *num_vars == 0 {
            fields |= Needs::GTS;
        }
        // The chromosome and the position are one field, and a block that
        // holds one of the two columns and not the other is a block that
        // `check` does not look at: it holds neither field.
        if chrom.is_some() && pos.is_some() {
            fields |= Needs::CHROM_POS;
        }
        if id.is_some() {
            fields |= Needs::ID;
        }
        if alleles.is_some() {
            fields |= Needs::ALLELES;
        }
        if qual.is_some() {
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
        let Block {
            num_vars: _,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = self;
        [
            !gts.is_empty(),
            chrom.is_some(),
            pos.is_some(),
            id.is_some(),
            alleles.is_some(),
            qual.is_some(),
        ]
    }

    /// How many alleles one variant of the block holds, the individuals
    /// times the ploidy, or the error of a block that cannot be. They are
    /// alleles and not genotypes: a genotype is the `ploidy` alleles of one
    /// individual, as `docs/glossary.md` has it.
    fn alleles_per_var(&self) -> Result<usize> {
        self.num_individuals
            .checked_mul(self.ploidy)
            .ok_or(Error::BlockTooLarge {
                num_vars_per_block: self.num_vars,
                num_individuals: self.num_individuals,
                ploidy: self.ploidy,
                // The block is here, so its size is one that a reader was
                // given and took: no machine allocated it, and what a
                // caller does about it is ask for fewer variants.
                size: BlockSize::AskedFor,
            })
    }

    /// The views of the variants of the block, in order.
    ///
    /// A view allocates nothing: its genotypes are a slice of `gts` and its
    /// other fields are read out of the columns.
    ///
    /// A block whose arrays are not of its size gives the views up to the
    /// first variant that is not in them and then stops, with no error, so
    /// a consumer whose blocks come from a reader with no [`Reblock`] and
    /// no binding crate in between calls [`Block::check`] first: both of
    /// those check every block they pass on.
    pub fn variants(&self) -> impl Iterator<Item = VariantRef<'_>> {
        (0..self.num_vars).map_while(|var| self.variant(var))
    }

    /// The view of the variant `var` of the block, counted from 0, or
    /// `None` when the block has no such variant.
    #[must_use]
    pub fn variant(&self, var: usize) -> Option<VariantRef<'_>> {
        let Block {
            num_vars,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = self;
        if var >= *num_vars {
            return None;
        }
        let gts = match gts.is_empty() {
            // A block built without the genotypes: every view has none.
            true => &[][..],
            false => {
                let alleles_per_var = self.alleles_per_var().ok()?;
                let start = var.checked_mul(alleles_per_var)?;
                let end = start.checked_add(alleles_per_var)?;
                gts.get(start..end)?
            }
        };
        Some(VariantRef::new(
            gts,
            chrom.as_ref().and_then(|column| column.get(var)).copied(),
            pos.as_ref().and_then(|column| column.get(var)).copied(),
            id.as_ref()
                .and_then(|column| column.get(var))
                .map(String::as_str),
            qual.as_ref().and_then(|column| column.get(var)).copied(),
            alleles.as_ref().map(|column| (column, var)),
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
        let alleles_per_var = self.alleles_per_var()?;
        let Block {
            num_vars,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = self;
        if !gts.is_empty() {
            let mut write = 0usize;
            for (var, keep_it) in keep.iter().enumerate() {
                if !*keep_it {
                    continue;
                }
                // Every one of these is a place in `gts`, which `check`
                // just said holds num_vars x alleles_per_var alleles.
                let read = var.saturating_mul(alleles_per_var);
                let end = read.saturating_add(alleles_per_var);
                if end <= gts.len() {
                    gts.copy_within(read..end, write);
                }
                write = write.saturating_add(alleles_per_var);
            }
            gts.truncate(write);
        }
        retain_in_column(chrom.as_mut(), keep);
        retain_in_column(pos.as_mut(), keep);
        retain_in_column(id.as_mut(), keep);
        retain_in_column(qual.as_mut(), keep);
        if let Some(alleles) = alleles.as_mut() {
            alleles.retain_vars(keep);
        }
        *num_vars = keep.iter().filter(|keep_it| **keep_it).count();
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
        let alleles_per_var = self.alleles_per_var()?;
        let Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = self;
        let alleles_of_the_block =
            num_vars
                .checked_mul(alleles_per_var)
                .ok_or(Error::BlockTooLarge {
                    num_vars_per_block: *num_vars,
                    num_individuals: *num_individuals,
                    ploidy: *ploidy,
                    size: BlockSize::AskedFor,
                })?;
        if !gts.is_empty() && gts.len() != alleles_of_the_block {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "gts",
                found: gts.len(),
                expected: alleles_of_the_block,
            });
        }
        let lengths = [
            ("chrom", chrom.as_ref().map(Vec::len)),
            ("pos", pos.as_ref().map(Vec::len)),
            ("id", id.as_ref().map(Vec::len)),
            ("qual", qual.as_ref().map(Vec::len)),
            ("alleles", alleles.as_ref().map(AllelesColumn::num_vars)),
        ];
        for (array, length) in lengths {
            if let Some(length) = length
                && length != *num_vars
            {
                return Err(Error::BlockArrayOfAnotherSize {
                    array,
                    found: length,
                    expected: *num_vars,
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
/// - Every block it gives passes [`Block::check`]: its genotypes are its
///   variants times its individuals times its ploidy, or none, and every
///   column it holds has one entry for each variant.
/// - Every block has the individuals of `individuals` and the ploidy of
///   `ploidy`, the same for every block of one pass.
/// - The chromosomes and the positions are two columns and one field: a
///   block holds both or neither, since [`Block::fields`] reports neither
///   when one of the two is missing.
///
/// [`Reblock`] gives each of those of the blocks of its source, and both
/// binding crates check a block before its genotypes cross to numpy or to
/// an `Int8Array`. A consumer that takes its blocks from a reader without
/// one of those in between calls [`Block::check`] itself before it walks
/// [`Block::variants`], which on a block whose arrays are too short stops
/// early and says nothing.
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

    /// The kind and the counts of every filter between this reader and its
    /// source, this one first when it is a filter.
    ///
    /// The kind is the name the counts have for a Python or a TypeScript
    /// user, `"missing_data"`, `"maf"` or `"obs_het"`. A source gives none,
    /// and a reader over another reader gives what its source gives, with
    /// its own before them when it is a filter. Whoever starts a pass keeps
    /// the chain of readers and reads the counts from it when the pass
    /// ends, as `docs/specs/filters.md` has it.
    ///
    /// The method has no default, so that a reader over another reader that
    /// forgets to pass on the counts of its source does not compile.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)>;
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

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        (**self).filtering_stats()
    }
}

/// A reader that is borrowed and not taken, which is how a consumer reads
/// a pass whose chain of readers its caller keeps.
///
/// `docs/specs/filters.md` has the counts of the filters of a pass read
/// from that chain when the consumer returns, so the consumer cannot own
/// it: `write_vars` of `docs/specs/io_vars.md` takes `&mut reader` and the
/// binding crate that built the chain reads the counts from it afterwards.
impl<R: BlockReader + ?Sized> BlockReader for &mut R {
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

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        (**self).filtering_stats()
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
/// ones that were waiting. Cutting copies out the rows that leave, into a
/// block allocated for them, and the rest stays in the block that waits,
/// with the row it starts at: so a block of 10000 variants cut into blocks
/// of 100 copies each row once and not once for every cut before it, and
/// the block that is given holds the memory of its own rows and not of the
/// block it was cut from.
pub struct Reblock<R: BlockReader> {
    reader: R,
    num_vars_per_block: usize,
    /// Whether that size is the one its caller asked for, which the error
    /// of a block the machine has no memory for names.
    size: BlockSize,
    /// The individuals and the ploidy of the source, which every block it
    /// gives has to have for its rows to be joined with the others'.
    num_individuals: usize,
    ploidy: usize,
    /// The block whose variants have not all been given: the one that did
    /// not fill a block, or the one a cut is taking blocks out of. It holds
    /// one variant at least when it is there.
    waiting: Option<Block>,
    /// How many rows of the block that waits were given already. The rows
    /// before it are not read again and their ids were taken out of it.
    given: usize,
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
        let (num_vars_per_block, size) =
            size_of_the_blocks(num_vars_per_block, num_individuals, ploidy)?;
        Ok(Reblock {
            reader,
            num_vars_per_block,
            size,
            num_individuals,
            ploidy,
            waiting: None,
            given: 0,
            finished: false,
        })
    }

    /// The error of a block whose memory the machine does not give.
    fn too_large(&self) -> Error {
        Error::BlockTooLarge {
            num_vars_per_block: self.num_vars_per_block,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
            size: self.size,
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
        let Block {
            num_vars: arrived,
            num_individuals: _,
            ploidy: _,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = block;
        let num_vars = waiting
            .num_vars
            .checked_add(arrived)
            .ok_or_else(|| self.too_large())?;
        try_extend(&mut waiting.gts, gts).map_err(|_| self.too_large())?;
        try_extend_column(waiting.chrom.as_mut(), chrom).map_err(|_| self.too_large())?;
        try_extend_column(waiting.pos.as_mut(), pos).map_err(|_| self.too_large())?;
        try_extend_column(waiting.id.as_mut(), id).map_err(|_| self.too_large())?;
        try_extend_column(waiting.qual.as_mut(), qual).map_err(|_| self.too_large())?;
        if let (Some(waiting), Some(arrived)) = (waiting.alleles.as_mut(), alleles.as_ref()) {
            waiting.try_append(arrived).map_err(|_| self.too_large())?;
        }
        waiting.num_vars = num_vars;
        Ok(())
    }

    /// How many rows of the block that waits have not been given yet.
    fn rows_waiting(&self) -> usize {
        self.waiting
            .as_ref()
            .map_or(0, |block| block.num_vars.saturating_sub(self.given))
    }

    /// The next `num_vars_per_block` rows of the block that waits, copied
    /// out of it, and the rest left where they are. A block that has the
    /// size and has given nothing goes through as it is, with no copy.
    fn cut(&mut self) -> Result<Option<Block>> {
        let Some(mut waiting) = self.waiting.take() else {
            return Ok(None);
        };
        if self.given == 0 && waiting.num_vars == self.num_vars_per_block {
            return Ok(Some(waiting));
        }
        let block = match take_rows(&mut waiting, self.given, self.num_vars_per_block) {
            Ok(block) => block,
            Err(_) => {
                self.finished = true;
                return Err(self.too_large());
            }
        };
        self.given = self.given.saturating_add(self.num_vars_per_block);
        if self.given < waiting.num_vars {
            self.waiting = Some(waiting);
        } else {
            self.given = 0;
        }
        Ok(Some(block))
    }

    /// What is left of the block that waits, as a block of its own: the
    /// block itself when no cut took rows out of it, and a copy of the rows
    /// that are left when one did.
    fn rest(&mut self) -> Result<Option<Block>> {
        let Some(mut waiting) = self.waiting.take() else {
            return Ok(None);
        };
        let left = waiting.num_vars.saturating_sub(self.given);
        if left == 0 {
            self.given = 0;
            return Ok(None);
        }
        let rest = match self.given {
            0 => waiting,
            given => match take_rows(&mut waiting, given, left) {
                Ok(rest) => rest,
                Err(_) => {
                    self.finished = true;
                    return Err(self.too_large());
                }
            },
        };
        self.given = 0;
        Ok(Some(rest))
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
            if self.rows_waiting() >= self.num_vars_per_block {
                return self.cut();
            }
            let block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    // What is left of the block that waits is the last
                    // block, shorter than the size that was asked for.
                    return self.rest();
                }
                Err(error) => {
                    // The source ends at its error and is not called
                    // again, and what was waiting is lost with it.
                    self.finished = true;
                    self.waiting = None;
                    self.given = 0;
                    return Err(error);
                }
            };
            if let Err(error) = self.taken(&block) {
                self.finished = true;
                self.waiting = None;
                self.given = 0;
                return Err(error);
            }
            let same_columns = self
                .waiting
                .as_ref()
                .is_some_and(|waiting| waiting.columns() == block.columns());
            if self.waiting.is_none() {
                self.waiting = Some(block);
                self.given = 0;
                continue;
            }
            if !same_columns {
                // The columns of the source changed, which a change of
                // `Needs` in the middle of a pass does: what is left of the
                // block that waits is given as a shorter block.
                let rest = self.rest().inspect_err(|_| {
                    self.finished = true;
                })?;
                self.waiting = Some(block);
                self.given = 0;
                return Ok(rest);
            }
            match self.rest() {
                Err(error) => {
                    self.finished = true;
                    return Err(error);
                }
                Ok(None) => {
                    self.waiting = Some(block);
                    self.given = 0;
                }
                Ok(Some(mut waiting)) => {
                    if let Err(error) = self.join(&mut waiting, block) {
                        self.finished = true;
                        return Err(error);
                    }
                    self.waiting = Some(waiting);
                    self.given = 0;
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

    /// The counts of the source: `reblock` takes no variant out, so it adds
    /// none of its own.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
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

/// `count` rows of `block` from the row `from` on, in a block of their own
/// that holds the columns and the individuals of the one they came from and
/// the memory of those rows alone.
///
/// `block` keeps its rows and its size: the caller of a cut says with the
/// row it starts at which of them it has given away. The ids are the one
/// column that is moved and not copied, an empty text left in the row that
/// leaves, so that a cut allocates nothing for each row.
///
/// # Errors
///
/// When the machine does not give the memory of the new block.
fn take_rows(
    block: &mut Block,
    from: usize,
    count: usize,
) -> std::result::Result<Block, TryReserveError> {
    // The block passed `check` when it was taken, so its arrays hold its
    // variants and these places are in them.
    let alleles_per_var = block.num_individuals.saturating_mul(block.ploidy);
    let mut gts = Vec::new();
    if !block.gts.is_empty() {
        let start = from.saturating_mul(alleles_per_var);
        let end = start.saturating_add(count.saturating_mul(alleles_per_var));
        let rows = block.gts.get(start..end).unwrap_or(&[]);
        gts.try_reserve_exact(rows.len())?;
        gts.extend_from_slice(rows);
    }
    let chrom = copied_rows(block.chrom.as_ref(), from, count)?;
    let pos = copied_rows(block.pos.as_ref(), from, count)?;
    let qual = copied_rows(block.qual.as_ref(), from, count)?;
    let alleles = match block.alleles.as_ref() {
        Some(column) => Some(column.try_rows(from, count)?),
        None => None,
    };
    let id = taken_rows(block.id.as_mut(), from, count)?;
    Ok(Block {
        num_vars: count,
        num_individuals: block.num_individuals,
        ploidy: block.ploidy,
        gts,
        chrom,
        pos,
        id,
        alleles,
        qual,
    })
}

/// `count` entries of `column` from `from` on, copied into a column of
/// their own whose memory is asked for with `try_reserve`: `Vec::to_vec`
/// allocates without asking, and a size that a caller wrote reaches it.
fn copied_rows<T: Copy>(
    column: Option<&Vec<T>>,
    from: usize,
    count: usize,
) -> std::result::Result<Option<Vec<T>>, TryReserveError> {
    let Some(column) = column else {
        return Ok(None);
    };
    let from = from.min(column.len());
    let end = from.saturating_add(count).min(column.len());
    let entries = column.get(from..end).unwrap_or(&[]);
    let mut rows = Vec::new();
    rows.try_reserve_exact(entries.len())?;
    rows.extend_from_slice(entries);
    Ok(Some(rows))
}

/// `count` ids of `column` from `from` on, moved into a column of their
/// own, each leaving an empty text where it was: an id is a `String`, and
/// copying one is an allocation for one row.
fn taken_rows(
    column: Option<&mut Vec<String>>,
    from: usize,
    count: usize,
) -> std::result::Result<Option<Vec<String>>, TryReserveError> {
    let Some(column) = column else {
        return Ok(None);
    };
    let from = from.min(column.len());
    let end = from.saturating_add(count).min(column.len());
    let mut rows = Vec::new();
    rows.try_reserve_exact(end.saturating_sub(from))?;
    for index in from..end {
        match column.get_mut(index) {
            Some(id) => rows.push(std::mem::take(id)),
            None => break,
        }
    }
    Ok(Some(rows))
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{
        AllelesColumn, Block, BlockReader, BlockSize, FIELD_NAMES, FIELDS_OF_THE_NAMES,
        GENOTYPES_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK, MIN_NUM_VARS_PER_BLOCK, Reblock,
        check_the_size_of_a_block, default_num_vars_per_block, needs_of_the_fields,
        size_of_the_blocks,
    };
    use crate::error::{Error, Result};
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs, VariantRef};

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
            num_vars_per_block: None,
        }
    }

    /// The same options with the blocks of the size that a test asks for.
    fn in_blocks_of(options: VcfOptions, num_vars_per_block: Option<usize>) -> VcfOptions {
        VcfOptions {
            num_vars_per_block,
            ..options
        }
    }

    /// A reader over one of the reference VCFs, which gives its blocks of
    /// `num_vars_per_block` variants holding `needs` and the genotypes.
    fn reader_over(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<BufReader<File>> {
        let options = in_blocks_of(options, num_vars_per_block);
        let mut reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        reader.set_needs(needs.union(Needs::GTS));
        reader
    }

    /// Every block a reader gives, until it has no more or it fails.
    fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
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
        let mut reader = reader_over(name, options, needs, num_vars_per_block);
        match blocks_of(&mut reader) {
            Ok(blocks) => blocks,
            Err(error) => panic!("{name}: the reader stopped at {error}"),
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

    /// The genotypes, the chromosome and the position of one variant,
    /// which is what the blocks, joined, have to hold whatever their
    /// size.
    #[derive(Debug, PartialEq, Eq)]
    struct Site {
        chrom: u32,
        pos: u64,
        gts: Vec<i8>,
    }

    /// The variants of a VCF in one block, which is what the blocks of
    /// every other size are compared with.
    fn sites_read_in_one_block(name: &str, options: VcfOptions) -> Vec<Site> {
        let blocks = blocks_read(
            name,
            options,
            Needs::CHROM_POS,
            Some(MAX_NUM_VARS_PER_BLOCK),
        );
        assert_eq!(blocks.len(), 1, "{name} is not one block");
        sites_of(&blocks)
    }

    /// One variant with every field it has, which is what a block holds
    /// when every field was asked for. The quality is `None` for a variant
    /// that has none, where the column of a block has a NaN, so that the
    /// comparison of two variants compares no float with another.
    #[derive(Debug, PartialEq)]
    struct FullVariant {
        chrom: u32,
        pos: u64,
        id: String,
        alleles: Vec<String>,
        qual: Option<f32>,
        gts: Vec<i8>,
    }

    /// Every field of every variant of a VCF in one block: what the blocks
    /// of that file have to hold, joined, whatever their size.
    fn variants_read_in_one_block(name: &str, options: VcfOptions) -> Vec<FullVariant> {
        let blocks = blocks_read(name, options, Needs::ALL, Some(MAX_NUM_VARS_PER_BLOCK));
        assert_eq!(blocks.len(), 1, "{name} is not one block");
        variants_of(&blocks)
    }

    /// Every field of every variant of the blocks, joined, read through
    /// their views.
    fn variants_of(blocks: &[Block]) -> Vec<FullVariant> {
        let mut variants = Vec::new();
        for block in blocks {
            block.check().expect("the block is of its size");
            for view in block.variants() {
                let num_alleles = view.num_alleles().expect("the alleles were asked for");
                variants.push(FullVariant {
                    chrom: view.chrom().expect("the chromosome"),
                    pos: view.pos().expect("the position"),
                    id: view.id().expect("the id").to_string(),
                    alleles: (0..num_alleles)
                        .map(|allele| view.allele(allele).unwrap_or("").to_string())
                        .collect(),
                    qual: view.qual().filter(|qual| !qual.is_nan()),
                    gts: view.gts().to_vec(),
                });
            }
        }
        variants
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
        let expected = sites_read_in_one_block("many.vcf", VcfOptions::default());
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
    fn a_reader_asked_for_the_genotypes_alone_gives_a_block_with_no_other_column() {
        // The VCF reader parses the chromosome and the position of every
        // line it gives a row to, and the block still has no column for
        // them.
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

    /// The error of a block the machine cannot give the memory for says
    /// which of the two sizes of a block it is about, because what a caller
    /// does differs: a size they asked for is one they lower, and the size
    /// popnei chose is one they replace with a size of their own.
    ///
    /// A `usize` is 64 bits natively and 32 in wasm, so the sizes that
    /// overflow it here are not the ones a user of the browser meets: what
    /// this checks is the three steps, which are the same on both.
    #[test]
    fn the_error_of_a_block_that_does_not_fit_says_whose_size_it_is() {
        // Two individuals of the ploidy 2 are four alleles in every
        // variant, so every size above a quarter of `usize::MAX` is one
        // that does not fit.
        let error = match size_of_the_blocks(Some(usize::MAX), 2, 2) {
            Ok(size) => panic!("the size {size:?} was taken"),
            Err(error) => error,
        };
        let message = error.to_string();
        assert!(matches!(error, Error::BlockTooLarge { .. }), "{message}");
        assert!(message.contains("ask for fewer variants"), "{message}");

        // A source of so many individuals that the size popnei chooses for
        // them, the smallest it chooses, does not fit either.
        let error = match size_of_the_blocks(None, usize::MAX / 4, 2) {
            Ok(size) => panic!("the size {size:?} was taken"),
            Err(error) => error,
        };
        let message = error.to_string();
        let Error::BlockTooLarge {
            num_vars_per_block, ..
        } = error
        else {
            panic!("the error is {message}");
        };
        assert_eq!(num_vars_per_block, MIN_NUM_VARS_PER_BLOCK);
        assert!(message.contains("popnei chose"), "{message}");
        assert!(message.contains("num_vars_per_block"), "{message}");
        assert!(!message.contains("ask for fewer variants"), "{message}");

        // And no block at all, which is neither of the two.
        let error = match size_of_the_blocks(Some(0), 2, 2) {
            Ok(size) => panic!("the size {size:?} was taken"),
            Err(error) => error,
        };
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");

        // A size that fits gives the number of the alleles of a full block,
        // which is what a reader reserves.
        let (num_vars_per_block, size) = size_of_the_blocks(Some(7), 2, 2).expect("the size");
        assert_eq!((num_vars_per_block, size), (7, BlockSize::AskedFor));
        assert_eq!(
            check_the_size_of_a_block(7, 2, 2, size).expect("the alleles"),
            28
        );
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

    /// A reader over a VCF written in a test.
    ///
    /// What a reader with no variant, and what an error in the middle of a
    /// block, give is tested where the reader is, in `io::vcf`.
    fn reader_over_text(
        vcf: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<Cursor<Vec<u8>>> {
        let options = in_blocks_of(options, num_vars_per_block);
        let mut reader = match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        };
        reader.set_needs(needs.union(Needs::GTS));
        reader
    }

    /// The columns of a block are asked of the machine when the block is
    /// started, and a size that a caller wrote reaches neither an abort nor
    /// a panic: the positions of a variant are 8 bytes whatever the
    /// individuals are, so a block of that many variants is refused for
    /// the memory of its columns although its genotypes were not asked
    /// for. It is made over the VCF reader, which is the reader that
    /// allocates blocks.
    #[test]
    fn a_block_of_more_memory_than_the_machine_gives_is_refused_when_it_is_started() {
        // Three individuals of the ploidy 2 are six alleles in a variant,
        // so a size of an eighth of `usize::MAX` is a multiplication that
        // does not carry over and a column that no machine gives: the
        // positions of those variants alone are eight times more bytes.
        let size = usize::MAX / 8;
        let vcf = vcf_of(&["chr1 100 rs1 A T . PASS . GT 0/0 0/1 1/1"]);
        let mut reader =
            reader_over_text(&vcf, VcfOptions::default(), Needs::CHROM_POS, Some(size));
        // The genotypes are not among what this block would hold, and the
        // columns of the chromosomes and the positions are.
        reader.set_needs(Needs::CHROM_POS);

        let error = match reader.next_block() {
            Ok(block) => panic!("the reader gave {block:?}"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((num_vars_per_block, num_individuals, ploidy), (size, 3, 2));

        assert!(reader.next_block().expect("no block").is_none());
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
            num_vars_per_block: None,
        };
        let mut reader = reader_over_text(&vcf, options, Needs::GTS, Some(2));
        let blocks = blocks_of(&mut reader).expect("the blocks");

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

    /// The two variants of the table of `differences.vcf` of that spec,
    /// three individuals of the ploidy 2 as well, whose alleles are of one
    /// byte and of several: a column that moved a text by a number of
    /// alleles where it should move it by a number of bytes gives the same
    /// answer for every allele of one byte, and not for these.
    const OF_SEVERAL_BYTES: [Row; 2] = [
        Row {
            pos: 50,
            id: "ms1",
            alleles: &["GTC", "G", "GTCT"],
            qual: Some(50.0),
            gts: [0, 1, 0, 2, MISSING, MISSING],
        },
        Row {
            pos: 60,
            id: "",
            alleles: &["A", "<DEL>", "*"],
            qual: None,
            gts: [0, 1, 2, 2, 0, 0],
        },
    ];

    /// A block built by hand from the rows of that table, with every
    /// column: this is what a reader of `cases.vcf` gives, and building it
    /// here is what lets the tests of the block and of `reblock` stand on
    /// the literals of the spec without a file.
    fn cases_block(rows: &[usize]) -> Block {
        let rows: Vec<&Row> = rows.iter().map(|row| &CASES[*row]).collect();
        block_of(&rows)
    }

    /// A block of the rows of either table, in the order they are given.
    fn block_of(rows: &[&Row]) -> Block {
        let mut gts = Vec::new();
        let mut chrom = Vec::new();
        let mut pos = Vec::new();
        let mut id = Vec::new();
        let mut qual = Vec::new();
        let mut alleles = AllelesColumn::with_num_vars(rows.len()).expect("the alleles");
        for row in rows {
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
    /// field of it.
    fn assert_view_is_the_row(view: &VariantRef<'_>, row: usize) {
        assert_view_is(view, &CASES[row], &format!("the row {row}"));
    }

    /// That a view holds every field of `expected`. The quality is compared
    /// as an `Option`, with the NaN of a variant that has none as `None`.
    fn assert_view_is(view: &VariantRef<'_>, expected: &Row, name: &str) {
        assert_eq!(view.gts(), expected.gts, "the genotypes of {name}");
        assert_eq!(view.chrom(), Some(0), "the chromosome of {name}");
        assert_eq!(view.pos(), Some(expected.pos), "the position of {name}");
        assert_eq!(view.id(), Some(expected.id), "the id of {name}");
        let qual = view.qual().filter(|qual| !qual.is_nan());
        assert_eq!(qual, expected.qual, "the quality of {name}");
        assert_eq!(
            alleles_of_view(view),
            expected.alleles,
            "the alleles of {name}"
        );
    }

    /// That the views of `blocks`, joined, are `rows`, in their order.
    fn assert_views_are_the_rows(blocks: &[Block], rows: &[&Row]) {
        let views: Vec<VariantRef<'_>> = blocks.iter().flat_map(Block::variants).collect();
        assert_eq!(views.len(), rows.len(), "how many variants the blocks hold");
        for (index, (view, row)) in views.iter().zip(rows).enumerate() {
            assert_view_is(view, row, &format!("the variant {index}, counted from 0"));
        }
        for block in blocks {
            block.check().expect("the block is of its size");
        }
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
        // The last variant of the block and the one after it, on a block
        // whose views have no genotypes: there the number of the variant
        // is the only bound there is.
        assert!(block.variant(3).is_some());
        assert!(block.variant(4).is_none());
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

        // Each of the five columns, one entry short of the variants of the
        // block: the error names the column that is short.
        for array in FIELD_NAMES {
            let mut block = cases_block(&[0, 1, 2, 3]);
            match array {
                "chrom" => drop(block.chrom.as_mut().expect("the chromosomes").pop()),
                "pos" => drop(block.pos.as_mut().expect("the positions").pop()),
                "id" => drop(block.id.as_mut().expect("the ids").pop()),
                // The column of the alleles counts its variants itself, so
                // one of three variants takes the place of one of four.
                "alleles" => block.alleles = cases_block(&[0, 1, 2]).alleles,
                _ => drop(block.qual.as_mut().expect("the qualities").pop()),
            }
            let error = match block.check() {
                Ok(()) => panic!("a block whose `{array}` is short passed the check"),
                Err(error) => error,
            };
            let Error::BlockArrayOfAnotherSize {
                array: named,
                found,
                expected,
            } = error
            else {
                panic!("the error of a short `{array}` is {error}");
            };
            assert_eq!((named, found, expected), (array, 3, 4));
        }
    }

    /// The arrays of a block are counted in variants, and a block of more
    /// variants than a `usize` holds genotypes for is the error of a block
    /// too large and not an overflow.
    #[test]
    fn a_block_of_more_genotypes_than_a_usize_holds_fails_check() {
        let of_the_size = |num_vars: usize, num_individuals: usize| Block {
            num_vars,
            num_individuals,
            ploidy: 2,
            gts: Vec::new(),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };

        // The individuals times the ploidy, the alleles of one variant.
        let error = match of_the_size(1, usize::MAX).check() {
            Ok(()) => panic!("the block passed the check"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_individuals, ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(num_individuals, usize::MAX);

        // And that times the variants of the block.
        let error = match of_the_size(usize::MAX, 2).check() {
            Ok(()) => panic!("the block passed the check"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block, ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(num_vars_per_block, usize::MAX);
    }

    /// `retain_vars` moves the rows by their place in the arrays, so it
    /// checks the block before it moves one: a block whose arrays are not
    /// of its size would be compacted with every row after the fault at
    /// the place of another.
    #[test]
    fn retain_vars_refuses_a_block_whose_arrays_are_not_of_its_size() {
        let mut block = cases_block(&[0, 1, 2, 3]);
        block.gts.pop();
        let error = match block.retain_vars(&[true, false, true, true]) {
            Ok(()) => panic!("the block was compacted"),
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
        assert_eq!((array, found, expected), ("gts", 23, 24));

        // The block is left as it was: its four variants and the allele it
        // was short.
        assert_eq!(block.num_vars, 4);
        assert_eq!(block.gts.len(), 23);
        assert_eq!(block.pos.as_ref().map(Vec::len), Some(4));
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
        // A block of no variants holds every field it was built with, the
        // genotypes among them, although its `gts` is empty.
        assert_eq!(block.fields(), Needs::ALL);
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
        // A user who gets this reads both counts.
        let message = error.to_string();
        assert!(message.contains('3'), "{message}");
        assert!(message.contains('4'), "{message}");

        // And five values for the same four variants, which is the other
        // way for a filter to be wrong.
        let error = match block.retain_vars(&[true; 5]) {
            Ok(()) => panic!("the block kept 5 values for 4 variants"),
            Err(error) => error,
        };
        let Error::KeepOfAnotherSize { found, num_vars } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((found, num_vars), (5, 4));

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
        /// What it reports of the filters between it and its source, so
        /// that a test sees whether a reader over it passes them on. A
        /// source has none, and this reader stands for a chain of filters
        /// when a test gives it some.
        filtering_stats: Vec<(&'static str, FilteringStats)>,
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
                filtering_stats: Vec::new(),
            }
        }

        /// The same reader, which reports `stats` as the counts of the
        /// filters between it and its source.
        fn reporting(
            blocks: Vec<Block>,
            stats: Vec<(&'static str, FilteringStats)>,
        ) -> GivenBlocks {
            GivenBlocks {
                filtering_stats: stats,
                ..GivenBlocks::of(blocks)
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

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            self.filtering_stats.clone()
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

    /// The counts of a chain of two filters, the outermost first, as they
    /// are on `many.vcf` with the missing data filter at 0.04 and the maf
    /// filter at 0.8 over it: 500 variants given and 215 kept, and then 215
    /// given and 163 kept, which `docs/specs/filters.md` has from bcftools
    /// 1.24 and pyNei. Here they are what a reader written for these tests
    /// reports, and no filter works them out.
    fn two_counts() -> Vec<(&'static str, FilteringStats)> {
        vec![
            (
                "maf",
                FilteringStats {
                    vars_processed: 215,
                    vars_kept: 163,
                },
            ),
            (
                "missing_data",
                FilteringStats {
                    vars_processed: 500,
                    vars_kept: 215,
                },
            ),
        ]
    }

    /// A source has no filter over it, and neither has a `reblock` over a
    /// source.
    #[test]
    fn a_source_and_a_reblock_over_it_give_no_filtering_stats() {
        let reader = GivenBlocks::of_one_variant_each();
        assert!(reader.filtering_stats().is_empty());
        let reblock = Reblock::new(reader, Some(3)).expect("the reblock");
        assert!(reblock.filtering_stats().is_empty());
    }

    /// `reblock` is a reader over a reader and passes on what its source
    /// reports, both counts and in the order it got them, before any block
    /// is taken and after the last one.
    #[test]
    fn reblock_gives_the_filtering_stats_of_its_source() {
        let mut reblock = Reblock::new(
            GivenBlocks::reporting(vec![cases_block(&[0, 1]), cases_block(&[2])], two_counts()),
            Some(2),
        )
        .expect("the reblock");
        assert_eq!(reblock.filtering_stats(), two_counts());
        let blocks = blocks_given(&mut reblock).expect("the blocks");
        assert_eq!(blocks.len(), 2);
        assert_eq!(reblock.filtering_stats(), two_counts());
    }

    /// A `Box<dyn BlockReader>`, which is how both binding crates hold the
    /// reader of a pass, gives what the reader inside it gives, and a
    /// `reblock` over a boxed reader gives the same.
    #[test]
    fn a_boxed_reader_gives_the_filtering_stats_of_the_reader_it_holds() {
        let reader: Box<dyn BlockReader> = Box::new(GivenBlocks::reporting(
            vec![cases_block(&[0, 1])],
            two_counts(),
        ));
        assert_eq!(reader.filtering_stats(), two_counts());
        let reblock = Reblock::new(reader, None).expect("the reblock");
        assert_eq!(reblock.filtering_stats(), two_counts());
        let boxed: Box<dyn BlockReader> = Box::new(GivenBlocks::of_one_variant_each());
        assert!(boxed.filtering_stats().is_empty());
    }

    /// A reader that is borrowed is the reader it borrows: it gives its
    /// blocks and its counts, what is asked of it reaches it, and it is
    /// still there when the borrow ends. That is how a consumer takes the
    /// chain of a pass whose caller keeps it, `write_vars` of
    /// `docs/specs/io_vars.md` over the chain a binding crate built.
    #[test]
    fn a_borrowed_reader_gives_the_blocks_and_the_counts_of_the_reader_it_borrows() {
        let mut reader =
            GivenBlocks::reporting(vec![cases_block(&[0, 1]), cases_block(&[2])], two_counts());
        {
            let mut reblock = Reblock::new(&mut reader, Some(2)).expect("the reblock");
            assert_eq!(reblock.filtering_stats(), two_counts());
            reblock.set_needs(Needs::GTS);
            let blocks = blocks_given(&mut reblock).expect("the blocks");
            assert_eq!(num_vars_of(&blocks), [2, 1]);
        }
        assert_eq!(reader.needs, Needs::GTS);
        assert_eq!(reader.filtering_stats(), two_counts());
        assert!(reader.left.is_empty());
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

    /// A cut copies out the rows it gives and leaves the rest where they
    /// are: the block that waits is the one the source gave, at the same
    /// address however many cuts were taken out of it, so a block of 500
    /// variants cut into blocks of 1 copies each row once and not once for
    /// every cut before it. And a block that was cut holds the memory of
    /// its own rows alone, which is what a Python user keeps when they hold
    /// its genotypes.
    #[test]
    fn a_cut_copies_the_rows_it_gives_and_leaves_the_rest_in_the_block_that_waits() {
        // The 500 variants of `many.vcf` in one block of the source.
        let source = source_over("many.vcf", every_variant(), Needs::CHROM_POS, Some(500));
        let mut reblock = Reblock::new(source, Some(1)).expect("the reblock");

        let mut blocks = Vec::new();
        let mut addresses_of_the_rest = Vec::new();
        while let Some(block) = reblock.next_block().expect("a block") {
            if let Some(waiting) = reblock.waiting.as_ref() {
                addresses_of_the_rest.push(waiting.gts.as_ptr().addr());
            }
            blocks.push(block);
        }

        assert_eq!(blocks.len(), 500);
        for block in &blocks {
            // 1 variant x 50 individuals x 2 alleles, and no more memory
            // than those.
            assert_eq!(block.gts.len(), 100);
            assert_eq!(block.gts.capacity(), 100);
            assert_eq!(block.pos.as_ref().map(Vec::capacity), Some(1));
        }
        // The 499 blocks that were given while rows were still waiting:
        // the same block of the source, neither copied nor allocated
        // again.
        assert_eq!(addresses_of_the_rest.len(), 499);
        let first = addresses_of_the_rest[0];
        assert!(
            addresses_of_the_rest
                .iter()
                .all(|address| *address == first),
            "the rest of the block was copied at a cut"
        );
    }

    /// The six rows of the two tables, the ones whose alleles are of
    /// several bytes among the ones whose alleles are of one: what a cut, a
    /// join and `retain_vars` have to carry. They are interleaved so that
    /// the alleles before any row are a different number of bytes and of
    /// alleles, which is what tells a shift in bytes from one in alleles.
    fn rows_of_both_tables() -> Vec<&'static Row> {
        vec![
            &CASES[0],
            &OF_SEVERAL_BYTES[0],
            &CASES[1],
            &OF_SEVERAL_BYTES[1],
            &CASES[2],
            &CASES[3],
        ]
    }

    /// A cut carries every column of the rows that leave, the alleles
    /// among them: a column that moved a text by a number of alleles and
    /// not by a number of bytes gives the alleles of one byte right and
    /// these wrong.
    #[test]
    fn a_cut_carries_every_column_of_the_rows_it_gives() {
        let rows = rows_of_both_tables();
        // One block of the source with the six rows, cut into blocks of 4:
        // the cut falls inside the rows whose alleles are of one byte, and
        // the block that is left holds the ones of several.
        let source = GivenBlocks::of(vec![block_of(&rows)]);
        let mut reblock = Reblock::new(source, Some(4)).expect("the reblock");
        let blocks = blocks_given(&mut reblock).expect("the blocks");

        assert_eq!(num_vars_of(&blocks), [4, 2]);
        assert_views_are_the_rows(&blocks, &rows);
    }

    /// A join carries every column of the block that arrives, the alleles
    /// among them, after the ones that were waiting.
    #[test]
    fn a_join_carries_every_column_of_the_blocks_it_joins() {
        let rows = rows_of_both_tables();
        // Six blocks of one variant each, joined into one of six.
        let blocks = rows.iter().map(|row| block_of(&[row])).collect();
        let mut reblock = Reblock::new(GivenBlocks::of(blocks), Some(6)).expect("the reblock");
        let given = blocks_given(&mut reblock).expect("the blocks");

        assert_eq!(num_vars_of(&given), [6]);
        assert_views_are_the_rows(&given, &rows);
    }

    /// `retain_vars` carries every column of the variants that stay, the
    /// alleles among them, and the variants it drops before them are what
    /// move the texts of the ones that stay.
    #[test]
    fn retain_vars_keeps_the_alleles_of_the_variants_that_stay() {
        let rows = rows_of_both_tables();
        let mut block = block_of(&rows);
        // A row whose alleles are of several bytes goes, and a row whose
        // alleles are of one, so the texts of the rows that stay move over
        // both.
        block
            .retain_vars(&[true, false, true, true, false, true])
            .expect("the variants to keep");

        assert_eq!(block.num_vars, 4);
        let kept: Vec<&Row> = vec![rows[0], rows[2], rows[3], rows[5]];
        assert_views_are_the_rows(std::slice::from_ref(&block), &kept);
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

    /// The chromosome and the position are one field and two columns, so
    /// `fields` says the same of a block that holds the chromosomes alone
    /// and of one that holds the positions alone: neither field. Their
    /// rows are still not rows that can be joined.
    #[test]
    fn reblock_does_not_join_blocks_whose_columns_differ_in_the_same_fields() {
        let mut with_the_chromosomes = cases_block(&[0]);
        with_the_chromosomes.pos = None;
        let mut with_the_positions = cases_block(&[1]);
        with_the_positions.chrom = None;
        assert_eq!(
            with_the_chromosomes.fields(),
            with_the_positions.fields(),
            "the two blocks hold the same fields"
        );

        let blocks = vec![with_the_chromosomes, with_the_positions];
        let mut reblock = Reblock::new(GivenBlocks::of(blocks), Some(2)).expect("the reblock");
        let given = blocks_given(&mut reblock).expect("the blocks");

        assert_eq!(num_vars_of(&given), [1, 1]);
        let first = given[0].variant(0).expect("the variant");
        assert_eq!((first.chrom(), first.pos()), (Some(0), None));
        let second = given[1].variant(0).expect("the variant");
        assert_eq!((second.chrom(), second.pos()), (None, Some(200)));
    }

    /// A block of more genotypes than this machine addresses is refused
    /// where the `Reblock` is built and before a block is read, on a
    /// machine of 32 bit addresses as on one of 64.
    #[test]
    fn a_reblock_of_more_genotypes_than_a_usize_holds_is_refused_when_it_is_built() {
        // The three individuals of the ploidy 2 of the reader are six
        // alleles in every variant, so every size above a sixth of
        // `usize::MAX` is one.
        let error = match Reblock::new(GivenBlocks::of(Vec::new()), Some(usize::MAX)) {
            Ok(reblock) => panic!("the reblock was built: {reblock:?}"),
            Err(error) => error,
        };
        let Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals,
            ploidy,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (num_vars_per_block, num_individuals, ploidy),
            (usize::MAX, 3, 2)
        );
    }

    /// Both binding crates hold their reader boxed, so what is asked of
    /// the box reaches the reader inside it.
    #[test]
    fn what_is_asked_of_a_boxed_reader_reaches_the_reader_inside_it() {
        let mut reader = Box::new(GivenBlocks::of(vec![cases_block(&[0])]));
        BlockReader::set_needs(&mut reader, Needs::GTS | Needs::QUAL);
        assert_eq!(reader.needs, Needs::GTS | Needs::QUAL);
        assert_eq!(BlockReader::individuals(&reader), ["ind1", "ind2", "ind3"]);
        assert_eq!(BlockReader::ploidy(&reader), 2);
        assert_eq!(BlockReader::chroms(&reader).name(0), Some("chr1"));
        let block = BlockReader::next_block(&mut reader)
            .expect("the block")
            .expect("a block");
        assert_eq!(block.num_vars, 1);
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

    /// The blocks of one of the reference VCFs, which `reblock` takes as
    /// its source.
    fn source_over(
        name: &str,
        options: VcfOptions,
        needs: Needs,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<BufReader<File>> {
        reader_over(name, options, needs, num_vars_per_block)
    }

    /// An error loses the block it happened in, and `reblock` loses with
    /// it the variants it was keeping for its next block. It is the case
    /// of "What a reader of the rules would not guess" of
    /// `docs/specs/block.md`, over a VCF written here.
    #[test]
    fn a_wrong_line_after_two_hundred_and_fifty_variants_leaves_twenty_eight_blocks_of_seven() {
        let mut lines: Vec<String> = (1..=250)
            .map(|variant| format!("chr1 {variant}00 rs{variant} A T . PASS . GT 0/0 0/1 1/1"))
            .collect();
        // The genotype of the 251st variant is of the ploidy 4 under a
        // reader of the ploidy 2.
        lines.push("chr1 25100 rs251 A T . PASS . GT 0/0 0/1 0/0/1/1".to_string());
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let vcf = vcf_of(&lines);

        let source = reader_over_text(&vcf, VcfOptions::default(), Needs::CHROM_POS, Some(100));
        let mut reblock = Reblock::new(source, Some(7)).expect("the reblock");

        let mut blocks = Vec::new();
        let error = loop {
            match reblock.next_block() {
                Ok(Some(block)) => blocks.push(block),
                Ok(None) => panic!("the reblock ended with no error"),
                Err(error) => break error,
            }
        };
        // The two blocks of 100 variants that the source gave whole are 28
        // blocks of 7, and the 4 variants that were left over are lost
        // with the error.
        assert_eq!(num_vars_of(&blocks), [7; 28]);
        assert!(
            matches!(error, Error::VcfGenotypePloidy { .. }),
            "the error is {error}"
        );
        assert!(reblock.next_block().expect("no block").is_none());
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
        // Every field of every variant of the file, read one variant at a
        // time by the VCF reader: what the blocks hold, joined, whatever
        // their size.
        let expected = variants_read_in_one_block("many.vcf", every_variant());
        assert_eq!(expected.len(), 500);
        for (num_vars_per_block, sizes_expected) in sizes {
            let source = source_over("many.vcf", every_variant(), Needs::ALL, Some(100));
            let mut reblock = Reblock::new(source, num_vars_per_block).expect("the reblock");
            let blocks = match blocks_given(&mut reblock) {
                Ok(blocks) => blocks,
                Err(error) => panic!("blocks of {num_vars_per_block:?}: {error}"),
            };
            assert_eq!(
                num_vars_of(&blocks),
                sizes_expected,
                "the variants of the blocks of {num_vars_per_block:?}"
            );
            let given = variants_of(&blocks);
            assert_eq!(
                given.len(),
                expected.len(),
                "blocks of {num_vars_per_block:?}: how many variants"
            );
            for (index, (given, expected)) in given.iter().zip(&expected).enumerate() {
                assert_eq!(
                    given, expected,
                    "blocks of {num_vars_per_block:?}: the variant {index}, counted from 0"
                );
            }
        }
    }
}
