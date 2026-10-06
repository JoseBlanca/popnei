//! The density of the variants along the chromosomes: how many variants
//! fall in each window of a given width, the windows of a chromosome laid
//! end to end from the position 1.
//!
//! Window k of a chromosome, counted from 0, holds the positions from k x
//! `window_size` + 1 to (k + 1) x `window_size`, both included. The last
//! window of a chromosome with a length ends at the length, and the last of
//! one with no length is the one that holds its last variant.
//! "The density of the variants along the chromosomes" of
//! `docs/specs/stats.md` has the design.

use std::collections::HashMap;
use std::fmt;
use std::num::NonZeroU64;

use super::{AfterABlock, SoFar, no_variant_in_the_pass, nothing_after_a_block};
use crate::block::{Block, BlockReader};
use crate::error::{Error, Result};
use crate::filters::FilteringStats;
use crate::phases::{Phase, timed};
use crate::variant::{ChromTable, Needs};

/// The most windows a density has, over all its chromosomes, which
/// `docs/specs/stats.md` decides: 10 million windows are 3.1 billion bases,
/// a human genome, in windows of 310 base pairs, and each is a row of the
/// frame of Python, about 100 bytes of it with the name of its chromosome.
pub const MAX_NUM_WINDOWS: usize = 10_000_000;

/// Where the lengths of the chromosomes of a density came from, which the
/// error of a variant past a length names so that a user knows which
/// length to look at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LengthsFrom {
    /// The lengths the caller gave, `chrom_lengths`, which replace those of
    /// the source for every chromosome.
    ChromLengths,
    /// The lengths of the header of the source: the `##contig` lines of a
    /// VCF that have a `length`, and what a vars file keeps of them.
    Source,
}

impl fmt::Display for LengthsFrom {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LengthsFrom::ChromLengths => formatter.write_str("`chrom_lengths`"),
            LengthsFrom::Source => formatter.write_str("the header of the source"),
        }
    }
}

/// How many variants of a pass fall in each window of each chromosome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VarDensity {
    window_size: NonZeroU64,
    chroms: Vec<DensityOfChrom>,
    num_vars: u64,
    num_windows: usize,
}

/// The windows of one chromosome.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DensityOfChrom {
    /// The name of the chromosome, as the source writes it.
    pub name: String,
    /// The length of the chromosome, its last position counted from 1.
    /// `None` when neither the caller nor the source gave one.
    pub length: Option<u64>,
    /// How many variants fall in each window, from the first, which starts
    /// at the position 1.
    pub counts: Vec<u32>,
}

/// One window of a density, which is one row of the frame of Python.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DensityWindow<'a> {
    /// The name of the chromosome of the window.
    pub chrom: &'a str,
    /// The first position of the window, counted from 1.
    pub start: u64,
    /// The last position of the window, included.
    pub end: u64,
    /// How many variants of the pass are at a position from `start` to
    /// `end`.
    pub num_vars: u32,
}

impl VarDensity {
    /// The width of a window in base pairs, 1 or more.
    #[must_use]
    pub fn window_size(&self) -> u64 {
        self.window_size.get()
    }

    /// The chromosomes in the order of the result: those with a length in
    /// the order of the lengths, then those with variants and no length in
    /// the order their first variant came.
    #[must_use]
    pub fn chroms(&self) -> &[DensityOfChrom] {
        &self.chroms
    }

    /// How many variants the pass gave, which is the sum of every count.
    #[must_use]
    pub fn num_vars(&self) -> u64 {
        self.num_vars
    }

    /// How many windows the density has over all its chromosomes,
    /// [`MAX_NUM_WINDOWS`] at most.
    #[must_use]
    pub fn num_windows(&self) -> usize {
        self.num_windows
    }

    /// Every window with its chromosome, its start, its end and its count:
    /// the chromosomes in the order of [`VarDensity::chroms`], and the
    /// windows of each in the order of their positions.
    ///
    /// Window k of a chromosome starts at k x `window_size` + 1 and ends
    /// `window_size` - 1 positions later, except for the last window of a
    /// chromosome with a length, which ends at the length, and the last of
    /// one with no length whose full width would pass
    /// 18446744073709551615, the largest position a source holds, which
    /// ends there.
    pub fn windows(&self) -> impl Iterator<Item = DensityWindow<'_>> {
        let window_size = self.window_size;
        self.chroms
            .iter()
            .flat_map(move |chrom| windows_of(chrom, window_size))
    }
}

/// The windows of one chromosome, from the first.
///
/// Each window starts where the one before it ended, plus 1. None of the
/// additions passes the largest `u64` but the end of the last window of a
/// chromosome with no length, which then ends at `u64::MAX`, as
/// [`VarDensity::windows`] says: a window is made only for the positions
/// up to the length, or up to a variant, so the start of every window is a
/// position a source holds.
fn windows_of(
    chrom: &DensityOfChrom,
    window_size: NonZeroU64,
) -> impl Iterator<Item = DensityWindow<'_>> {
    // The width is 1 at least, so this takes nothing away from it but the
    // first position.
    let after_the_start = window_size.get().saturating_sub(1);
    chrom.counts.iter().scan(1_u64, move |start, num_vars| {
        let full_width = start.saturating_add(after_the_start);
        let end = chrom
            .length
            .map_or(full_width, |length| full_width.min(length));
        let window = DensityWindow {
            chrom: &chrom.name,
            start: *start,
            end,
            num_vars: *num_vars,
        };
        *start = end.saturating_add(1);
        Some(window)
    })
}

/// The density of the variants of one pass over `reader` in windows of
/// `window_size` base pairs, what `calc_var_density` gives in Python and
/// `calcVarDensity` in TypeScript.
///
/// The pass asks the reader for the chromosome and the position alone. The
/// lengths of the chromosomes are `chrom_lengths` when it is given, and
/// the lengths of `reader.header()` when it is `None`; lengths that are
/// given replace those of the source for every chromosome, so one they do
/// not name has no length. A chromosome with a length has its windows up to
/// the length, whether or not it has a variant; one with no length has them
/// up to the window of its last variant. A window with no variant has a
/// count of 0. The variants need not be sorted.
///
/// # Errors
///
/// - [`Error::VarDensityWindowSizeZero`] for a `window_size` of 0.
/// - [`Error::VarDensityChromLengthZero`] and
///   [`Error::VarDensityChromLengthTwice`] for a length of 0 and for a
///   chromosome named twice among the lengths.
/// - [`Error::VarDensityTooManyWindows`] for more than [`MAX_NUM_WINDOWS`]
///   windows, before the pass when the lengths give them and when the pass
///   reaches them otherwise.
/// - [`Error::VarDensityVarPastTheLength`] and
///   [`Error::VarDensityVarAtPositionZero`] for a variant that lies in no
///   window of its chromosome.
/// - [`Error::VarDensityWindowTooFull`] for a window of more than
///   `u32::MAX` variants.
/// - [`Error::FieldsNotInTheBlock`] for a block with no chromosome or no
///   position, which a source of no positions gives: a vars file written
///   without them, or variants built from an array of genotypes.
/// - [`Error::PassGaveNoVariant`] for a pass that gave no variant, whatever
///   the lengths.
/// - What [`Block::check`] refuses, [`Error::ReaderGaveABlockOfNoVariants`]
///   for a block of no variants and [`Error::VarDensityChromNameMissing`]
///   for a chromosome number the table of the reader has no name for, all
///   defects of a reader, and whatever the reader fails with.
pub fn calc_var_density<R: BlockReader + ?Sized>(
    reader: &mut R,
    window_size: u64,
    chrom_lengths: Option<&[(String, u64)]>,
) -> Result<VarDensity> {
    calc_var_density_with(
        reader,
        window_size,
        chrom_lengths,
        &mut nothing_after_a_block,
    )
}

/// [`calc_var_density`], with `after_a_block` called after each block the
/// pass adds, the last one too, with the density over the variants read so
/// far.
///
/// A chromosome with a length has all its windows from the first call; one
/// with no length has them up to the window of its last variant read so
/// far, so its windows grow from one call to the next, and a chromosome
/// with no length that no variant has reached yet is not in the result.
///
/// # Errors
///
/// Those of [`calc_var_density`], and the error `after_a_block` returns,
/// which ends the pass there.
pub fn calc_var_density_with<R: BlockReader + ?Sized>(
    reader: &mut R,
    window_size: u64,
    chrom_lengths: Option<&[(String, u64)]>,
    after_a_block: AfterABlock<'_, VarDensity>,
) -> Result<VarDensity> {
    let mut density = TheDensity::before_the_pass(&*reader, window_size, chrom_lengths)?;
    // The chromosome and the position are all a window needs, so a reader
    // over a file leaves the genotypes and every other column unparsed.
    reader.set_needs(Needs::CHROM_POS);
    // The blocks are read on this thread, without the read ahead of the
    // other passes of this module: the work on a block is one division and
    // one addition for each variant, and the pass has not been measured
    // yet, which the first measurement of "Speed" of `docs/specs/stats.md`
    // does.
    while let Some(block) = timed(Phase::NextBlock, || reader.next_block())? {
        timed(Phase::Work, || {
            density.count_the_block(&block, reader.chroms())
        })?;
        after_a_block(&DensitySoFar {
            density: &density,
            chain: &*reader,
        })?;
    }
    if density.num_vars == 0 {
        return Err(no_variant_in_the_pass(reader));
    }
    Ok(density.into_the_result())
}

/// What a pass of [`calc_var_density_with`] has counted after a block,
/// which builds the density from its counts when asked.
struct DensitySoFar<'pass, R: ?Sized> {
    density: &'pass TheDensity,
    /// The chain of the pass, which has just given the last block.
    chain: &'pass R,
}

impl<R: BlockReader + ?Sized> SoFar<VarDensity> for DensitySoFar<'_, R> {
    fn num_vars(&self) -> u64 {
        self.density.num_vars
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.chain.filtering_stats()
    }

    fn result(&self) -> Result<VarDensity> {
        Ok(self.density.the_result())
    }
}

/// The counts of a density while the pass goes, which a pass of
/// [`calc_var_density_with`] counts from block to block and
/// [`calc_variants_summary`](super::calc_variants_summary) counts too.
#[derive(Debug)]
pub(super) struct TheDensity {
    window_size: NonZeroU64,
    /// Where the lengths came from, which the error of a variant past one
    /// names.
    from: LengthsFrom,
    chroms: TheChroms,
    /// The windows of every chromosome of `chroms`, [`MAX_NUM_WINDOWS`] at
    /// most.
    num_windows: u64,
    pub(super) num_vars: u64,
}

/// The chromosomes of a density, with where each is found by its name and
/// by its number in the table of the reader.
#[derive(Debug)]
struct TheChroms {
    /// The chromosomes met so far, those with a length first.
    chroms: Vec<DensityOfChrom>,
    /// The place in `chroms` of each chromosome with a length, by its name.
    of_the_lengths: HashMap<String, usize>,
    /// The place in `chroms` of each chromosome number of the table of the
    /// reader, `None` for a number no variant of the pass has had yet.
    of_the_numbers: Vec<Option<usize>>,
}

impl TheDensity {
    /// The density of a pass over `reader` before any block, in windows of
    /// `window_size`, with the lengths of `chrom_lengths` when it is given
    /// and those of the header of `reader` when it is not.
    ///
    /// # Errors
    ///
    /// A `window_size` of 0, and those of the lengths that
    /// [`calc_var_density`] gives before the pass.
    pub(super) fn before_the_pass<R: BlockReader + ?Sized>(
        reader: &R,
        window_size: u64,
        chrom_lengths: Option<&[(String, u64)]>,
    ) -> Result<TheDensity> {
        let window_size = NonZeroU64::new(window_size).ok_or(Error::VarDensityWindowSizeZero)?;
        match chrom_lengths {
            Some(lengths) => {
                TheDensity::of_the_lengths(lengths, LengthsFrom::ChromLengths, window_size)
            }
            None => TheDensity::of_the_lengths(
                &reader.header().chrom_lengths,
                LengthsFrom::Source,
                window_size,
            ),
        }
    }

    /// The density before the pass: every chromosome of `lengths` with all
    /// its windows at 0, in the order of `lengths`.
    ///
    /// # Errors
    ///
    /// A length of 0, a chromosome named twice, and more windows than
    /// [`MAX_NUM_WINDOWS`] over all of them.
    fn of_the_lengths(
        lengths: &[(String, u64)],
        from: LengthsFrom,
        window_size: NonZeroU64,
    ) -> Result<TheDensity> {
        let mut chroms = TheChroms {
            chroms: Vec::with_capacity(lengths.len()),
            of_the_lengths: HashMap::with_capacity(lengths.len()),
            of_the_numbers: Vec::new(),
        };
        let mut num_windows: u64 = 0;
        for (name, length) in lengths {
            let Some(before_the_last) = length.checked_sub(1) else {
                return Err(Error::VarDensityChromLengthZero {
                    chrom: name.clone(),
                    from,
                });
            };
            if chroms
                .of_the_lengths
                .insert(name.clone(), chroms.chroms.len())
                .is_some()
            {
                return Err(Error::VarDensityChromLengthTwice {
                    chrom: name.clone(),
                    from,
                });
            }
            // The window of the last position, and one more for the window
            // 0: a length of u64::MAX in windows of 1 is u64::MAX windows,
            // which the addition reaches and does not pass.
            let of_the_chrom = (before_the_last / window_size).saturating_add(1);
            num_windows = num_windows.saturating_add(of_the_chrom);
            let Some(of_the_chrom) = below_the_bound(num_windows, of_the_chrom) else {
                return Err(too_many_windows(num_windows, window_size));
            };
            chroms.chroms.push(DensityOfChrom {
                name: name.clone(),
                length: Some(*length),
                counts: vec![0; of_the_chrom],
            });
        }
        Ok(TheDensity {
            window_size,
            from,
            chroms,
            num_windows,
            num_vars: 0,
        })
    }

    /// It counts every variant of `block` in the window of its position.
    ///
    /// # Errors
    ///
    /// Those of [`calc_var_density`] that a block gives.
    pub(super) fn count_the_block(&mut self, block: &Block, names: &ChromTable) -> Result<()> {
        // The chromosomes and the positions are walked together, so a
        // column of another length than the block would pair a variant with
        // the position of another.
        block.check()?;
        // Every reader of popnei gives one variant at least in a block and
        // no block when it has no more, so a block of none is a defect of
        // the reader, which the other two statistics refuse too.
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        let missing = Needs::CHROM_POS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        let (Some(numbers), Some(positions)) = (&block.chrom, &block.pos) else {
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::CHROM_POS,
            });
        };
        for (number, pos) in numbers.iter().zip(positions) {
            let chrom = self.chroms.chrom_of(*number, names)?;
            count_the_variant(
                chrom,
                *pos,
                self.window_size,
                self.from,
                &mut self.num_windows,
            )?;
        }
        // A pass of more than 18446744073709551615 variants reads more rows
        // than any source holds.
        self.num_vars = self.num_vars.saturating_add(super::num_vars_of(block));
        Ok(())
    }

    /// The density of the pass when it ends, which takes the counts.
    pub(super) fn into_the_result(self) -> VarDensity {
        let num_windows = self.num_windows_of_the_result();
        VarDensity {
            window_size: self.window_size,
            chroms: self.chroms.chroms,
            num_vars: self.num_vars,
            num_windows,
        }
    }

    /// The density over the variants counted so far, which copies the
    /// counts and leaves them to the pass.
    pub(super) fn the_result(&self) -> VarDensity {
        VarDensity {
            window_size: self.window_size,
            chroms: self.chroms.chroms.clone(),
            num_vars: self.num_vars,
            num_windows: self.num_windows_of_the_result(),
        }
    }

    fn num_windows_of_the_result(&self) -> usize {
        // At most MAX_NUM_WINDOWS, which every addition to it checked.
        usize::try_from(self.num_windows).unwrap_or(MAX_NUM_WINDOWS)
    }
}

impl TheChroms {
    /// The chromosome of the density that the chromosome `number` of the
    /// table of the reader is: the one of its length, or, the first time a
    /// variant of a chromosome with no length comes, a new one after all
    /// those the density has.
    ///
    /// # Errors
    ///
    /// When the table has no name for `number`.
    fn chrom_of(&mut self, number: u32, names: &ChromTable) -> Result<&mut DensityOfChrom> {
        let missing = Error::VarDensityChromNameMissing { number };
        let Ok(index) = usize::try_from(number) else {
            return Err(missing);
        };
        let place = if let Some(Some(place)) = self.of_the_numbers.get(index) {
            *place
        } else {
            let Some(name) = names.name(number) else {
                return Err(missing);
            };
            let place = if let Some(place) = self.of_the_lengths.get(name) {
                *place
            } else {
                let place = self.chroms.len();
                self.chroms.push(DensityOfChrom {
                    name: name.to_owned(),
                    length: None,
                    counts: Vec::new(),
                });
                place
            };
            // The table has a name for `number`, so it holds more names than
            // `index`: the vector grows to the table and no further, which is
            // memory the reader holds already.
            if self.of_the_numbers.len() < names.len() {
                self.of_the_numbers.resize(names.len(), None);
            }
            if let Some(entry) = self.of_the_numbers.get_mut(index) {
                *entry = Some(place);
            }
            place
        };
        // Every place of `of_the_numbers` and of `of_the_lengths` is that of
        // a chromosome pushed into `chroms`, so this finds it.
        self.chroms.get_mut(place).ok_or(missing)
    }
}

/// `added` as the number of windows of a vector, when the density, which
/// has `num_windows` with them, has [`MAX_NUM_WINDOWS`] at most.
///
/// A `usize` is 32 bits in wasm, so the numbers are held in a `u64` until
/// they are known to be below the bound, which fits in both.
fn below_the_bound(num_windows: u64, added: u64) -> Option<usize> {
    let fits = usize::try_from(num_windows).ok()?;
    if fits > MAX_NUM_WINDOWS {
        return None;
    }
    // At most `num_windows`, which fitted.
    usize::try_from(added).ok()
}

/// The error of a density of `num_windows` windows at least, which
/// saturates at `u64::MAX`: that is as many as it says or more.
fn too_many_windows(num_windows: u64, window_size: NonZeroU64) -> Error {
    Error::VarDensityTooManyWindows {
        num_windows,
        window_size: window_size.get(),
        largest: MAX_NUM_WINDOWS,
    }
}

/// It adds the variant at `pos` of `chrom` to the count of its window, with
/// the windows of a chromosome with no length made up to that window, and
/// `num_windows` of the density counted up with them.
///
/// # Errors
///
/// A position of 0 or past the length, more windows than
/// [`MAX_NUM_WINDOWS`], and a window of more than `u32::MAX` variants.
fn count_the_variant(
    chrom: &mut DensityOfChrom,
    pos: u64,
    window_size: NonZeroU64,
    from: LengthsFrom,
    num_windows: &mut u64,
) -> Result<()> {
    let Some(from_the_first) = pos.checked_sub(1) else {
        return Err(Error::VarDensityVarAtPositionZero {
            chrom: chrom.name.clone(),
        });
    };
    if let Some(length) = chrom.length
        && pos > length
    {
        return Err(Error::VarDensityVarPastTheLength {
            chrom: chrom.name.clone(),
            pos,
            length,
            from,
        });
    }
    let window = from_the_first / window_size;
    if let Ok(index) = usize::try_from(window)
        && let Some(count) = chrom.counts.get_mut(index)
    {
        let Some(with_it) = count.checked_add(1) else {
            return Err(window_too_full(chrom, window, window_size));
        };
        *count = with_it;
        return Ok(());
    }
    // A window past the vector is past the length of every chromosome that
    // has one, whose windows were all made before the pass, so it is of a
    // chromosome with no length. It is made here, with the ones between it
    // and the last one there was, and it holds this variant. `window` is at
    // most u64::MAX - 1, the last position in windows of 1, so one more is
    // a u64.
    let made = u64::try_from(chrom.counts.len()).unwrap_or(u64::MAX);
    let added = window.saturating_add(1).saturating_sub(made);
    let with_them = num_windows.saturating_add(added);
    if below_the_bound(with_them, added).is_none() {
        return Err(too_many_windows(with_them, window_size));
    }
    // Below the bound, as the windows of the density are.
    let Ok(index) = usize::try_from(window) else {
        return Err(too_many_windows(with_them, window_size));
    };
    *num_windows = with_them;
    chrom.counts.resize(index, 0);
    chrom.counts.push(1);
    Ok(())
}

/// The error of the window `window` of `chrom`, counted from 0, which holds
/// `u32::MAX` variants and was given one more.
fn window_too_full(chrom: &DensityOfChrom, window: u64, window_size: NonZeroU64) -> Error {
    // The window holds a variant, so its start is a position of the source
    // and none of these saturates but the end of a chromosome with no
    // length, which `VarDensity::windows` ends at u64::MAX as well.
    let start = window.saturating_mul(window_size.get()).saturating_add(1);
    let full_width = start.saturating_add(window_size.get().saturating_sub(1));
    Error::VarDensityWindowTooFull {
        chrom: chrom.name.clone(),
        start,
        end: chrom
            .length
            .map_or(full_width, |length| full_width.min(length)),
        largest: u32::MAX,
    }
}

#[cfg(test)]
mod tests;
