//! The filter by regions of `docs/specs/filters.md`: it keeps the variants
//! inside the regions of a BED file, or, with `exclude`, those outside all
//! of them.
//!
//! [`Regions`] reads the BED, plain or gzipped, and holds its regions
//! joined and sorted within each chromosome; [`RegionSelection`] is those
//! regions with which side of them is kept, and is what a source is offered
//! with [`BlockReader::skip_outside`]; [`RegionFilter`] keeps the variants
//! of a block that the selection keeps, with its two counts; and
//! [`RegionsReader`] puts that filter over a reader.
//!
//! A region of BED counts the bases from 0 and leaves its end out, and a
//! region of popnei counts them from 1, as a VCF does, and holds both ends:
//! the line `chr1 0 2000` is the region from 1 to 2000, and a line of start
//! `s` and end `e` the region from `s + 1` to `e`. A variant is inside when
//! its position, the POS of its VCF, is in a region of its chromosome.

use std::collections::HashMap;
use std::fmt;
use std::io::Read;
use std::sync::Arc;

use flate2::read::MultiGzDecoder;

use crate::block::{Block, BlockReader, SourceHeader};
use crate::error::{Error, Result};
use crate::filters::FilteringStats;
use crate::variant::{ChromTable, Needs};

/// The two bytes every gzipped file starts with.
const GZIP_BYTES: [u8; 2] = [0x1f, 0x8b];

/// The three bytes of the byte order mark of UTF-8, which are dropped from
/// the start of the text of a BED.
const UTF8_BYTE_ORDER_MARK: &[u8] = b"\xef\xbb\xbf";

/// The first words of the two lines the tools of the UCSC genome browser
/// write, which are not regions and which bedtools skips too. A comment and
/// an empty line are skipped as well.
const WORDS_OF_THE_LINES_SKIPPED: [&[u8]; 2] = [b"track", b"browser"];

/// The kind of the filter that keeps the variants inside the regions, and
/// of the one that keeps those outside them: the names its counts and its
/// step have for a Python and a TypeScript user.
const THE_KIND_INSIDE: &str = "regions";
const THE_KIND_OUTSIDE: &str = "excluded_regions";

/// One region, counted from 1 and with both ends in it, as a region of
/// popnei is: `first` is the start of its line of BED plus 1 and `last` is
/// the end of that line.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Region {
    first: u64,
    last: u64,
}

/// The regions of a BED file, joined where they overlap or touch and
/// sorted within each chromosome, which the filter by regions keeps the
/// variants inside or outside of.
///
/// They are read once, when a user adds the step, and shared behind an
/// `Arc` by the step and by every pass built from it.
#[derive(Clone, PartialEq, Eq)]
pub struct Regions {
    /// The regions of each chromosome, by its name as the BED writes it,
    /// sorted by their first position, none of them overlapping or touching
    /// the next.
    of_each_chrom: HashMap<Vec<u8>, Vec<Region>>,
    /// How many regions there are over every chromosome.
    num_regions: usize,
}

/// What is wrong with a line of a BED file, which [`Error::BedLine`] says
/// with the number of the line.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum BedLineProblem {
    /// The line has fewer than three columns separated by tabs: a region is
    /// the chromosome, the start and the end.
    FewerThanThreeColumns {
        /// How many columns separated by tabs it has.
        columns: usize,
        /// Whether it holds a space, which is what a line whose columns are
        /// separated by spaces holds, so that the message says BED
        /// separates them by tabs.
        with_spaces: bool,
    },
    /// The chromosome is empty: the line starts with a tab, and names no
    /// chromosome a variant can be on.
    EmptyChromosome,
    /// The start is not a whole number of 0 or more written in digits
    /// alone.
    StartNotAWholeNumber {
        /// The start as the line writes it.
        found: String,
    },
    /// The end is not a whole number of 0 or more written in digits
    /// alone.
    EndNotAWholeNumber {
        /// The end as the line writes it.
        found: String,
    },
    /// The start is a whole number above the largest of 64 bits.
    StartAboveTheLargest {
        /// The start as the line writes it.
        found: String,
    },
    /// The end is a whole number above the largest of 64 bits.
    EndAboveTheLargest {
        /// The end as the line writes it.
        found: String,
    },
    /// The start is not below the end. BED allows a start equal to its end
    /// for a point between two bases, which holds no position of popnei's.
    StartNotBelowEnd {
        /// The start of the line.
        start: u64,
        /// The end of the line.
        end: u64,
    },
}

impl fmt::Display for BedLineProblem {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BedLineProblem::FewerThanThreeColumns {
                columns,
                with_spaces,
            } => {
                let spaces = if *with_spaces {
                    ", and it holds spaces, which do not separate the columns of BED"
                } else {
                    ""
                };
                let columns_of = if *columns == 1 { "column" } else { "columns" };
                write!(
                    formatter,
                    "it has {columns} {columns_of} separated by tabs{spaces}; a region of BED is \
                     three columns at least, the chromosome, the start and the end, which BED \
                     separates by tabs"
                )
            }
            BedLineProblem::EmptyChromosome => write!(
                formatter,
                "its chromosome is empty, and a region is of a chromosome that has a name"
            ),
            BedLineProblem::StartNotAWholeNumber { found } => write!(
                formatter,
                "its start is `{found}`, and a start is a whole number of 0 or more, written in \
                 digits"
            ),
            BedLineProblem::EndNotAWholeNumber { found } => write!(
                formatter,
                "its end is `{found}`, and an end is a whole number of 0 or more, written in \
                 digits"
            ),
            BedLineProblem::StartAboveTheLargest { found } => write!(
                formatter,
                "its start, {found}, is above {largest}, the largest number of 64 bits",
                largest = u64::MAX
            ),
            BedLineProblem::EndAboveTheLargest { found } => write!(
                formatter,
                "its end, {found}, is above {largest}, the largest number of 64 bits",
                largest = u64::MAX
            ),
            BedLineProblem::StartNotBelowEnd { start, end } => write!(
                formatter,
                "its start, {start}, is not below its end, {end}; BED counts the bases from 0 \
                 and leaves the end out, so a region holds the positions from its start plus 1 \
                 to its end, and one whose start is not below its end holds none"
            ),
        }
    }
}

impl Regions {
    /// The regions of the BED text of `source`, plain or gzipped: gzipped
    /// when it starts with the two bytes of gzip, whatever the name of its
    /// file.
    ///
    /// A line that is empty, that starts with `#`, or whose first word is
    /// `track` or `browser`, is skipped. A byte order mark of UTF-8 at the
    /// start of the text is dropped. A carriage return before the end of a
    /// line is dropped,
    /// as bcftools and plink2 drop it. The columns after the third are not
    /// read. The lines are counted from 1 over the whole file, the skipped
    /// ones among them.
    ///
    /// # Errors
    ///
    /// [`Error::BedLine`], with the number of the line, for a line of fewer
    /// than three columns separated by tabs, an empty chromosome, a start or an end that is not
    /// a whole number of 0 or more that fits in 64 bits, and a start that
    /// is not below its end; [`Error::BedWithNoRegion`] for a BED with no
    /// line of a region; and [`Error::Io`] when `source` cannot be read or
    /// its gzip stream ends in the middle or is damaged.
    pub fn from_bed<S: Read>(mut source: S) -> Result<Regions> {
        let mut bytes = Vec::new();
        source.read_to_end(&mut bytes)?;
        Regions::from_bed_bytes(&bytes)
    }

    /// The regions of the BED whose bytes are `bytes`, plain or gzipped,
    /// which is [`Regions::from_bed`] for a caller that holds the bytes
    /// already: a plain BED is read where it is, with no copy of it, which
    /// in wasm, whose memory never shrinks, is memory the page keeps.
    ///
    /// # Errors
    ///
    /// Those of [`Regions::from_bed`], of which [`Error::Io`] only for a
    /// gzip stream that ends in the middle or is damaged.
    pub fn from_bed_bytes(bytes: &[u8]) -> Result<Regions> {
        if bytes.starts_with(&GZIP_BYTES) {
            // `MultiGzDecoder` goes on to the next member of the file when
            // one ends, so a BED that bgzip wrote, many members one after
            // another, is read to its end.
            let mut text = Vec::new();
            MultiGzDecoder::new(bytes).read_to_end(&mut text)?;
            return Regions::of_the_text(&text);
        }
        Regions::of_the_text(bytes)
    }

    /// The regions of the lines of `text`, which [`Regions::from_bed`]
    /// reads.
    fn of_the_text(text: &[u8]) -> Result<Regions> {
        // The byte order mark of UTF-8, which some editors of Windows write
        // at the start of a text, would be the start of the first name.
        let text = text.strip_prefix(UTF8_BYTE_ORDER_MARK).unwrap_or(text);
        let mut of_each_chrom: HashMap<Vec<u8>, Vec<Region>> = HashMap::new();
        let mut regions_read = false;
        for (index, line) in text.split(|byte| *byte == b'\n').enumerate() {
            // An index of a slice held in memory fits in a `u64` on every
            // target popnei builds for, so this does not saturate.
            let number = u64::try_from(index).unwrap_or(u64::MAX).saturating_add(1);
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if is_not_a_region(line) {
                continue;
            }
            let (chrom, region) = region_of_the_line(line).map_err(|problem| Error::BedLine {
                line: number,
                problem,
            })?;
            // The name is copied once for its chromosome and not for each of
            // its lines.
            match of_each_chrom.get_mut(chrom) {
                Some(regions) => regions.push(region),
                None => {
                    of_each_chrom.insert(chrom.to_vec(), vec![region]);
                }
            }
            regions_read = true;
        }
        if !regions_read {
            return Err(Error::BedWithNoRegion);
        }
        let mut num_regions: usize = 0;
        for regions in of_each_chrom.values_mut() {
            join(regions);
            num_regions = num_regions.saturating_add(regions.len());
        }
        Ok(Regions {
            of_each_chrom,
            num_regions,
        })
    }

    /// No region, which only the tests of the readers build: they offer it
    /// to [`BlockReader::skip_outside`] and look at the answer.
    #[cfg(test)]
    pub(crate) fn none_for_the_tests() -> Regions {
        Regions {
            of_each_chrom: HashMap::new(),
            num_regions: 0,
        }
    }

    /// How many regions are left once those that overlap or touch are
    /// joined, over every chromosome.
    #[must_use]
    pub fn num_regions(&self) -> usize {
        self.num_regions
    }

    /// Whether position `pos` of the chromosome named `chrom` is in a
    /// region: the name is compared with the one the BED writes byte for
    /// byte, and the position counts from 1, as the POS of a VCF does.
    #[must_use]
    pub fn contains(&self, chrom: &str, pos: u64) -> bool {
        self.of_the_chrom(chrom.as_bytes())
            .is_some_and(|regions| holds(regions, pos))
    }

    /// The regions of the chromosome of that name, sorted, or None when
    /// the BED names no region of it.
    fn of_the_chrom(&self, chrom: &[u8]) -> Option<&[Region]> {
        self.of_each_chrom.get(chrom).map(Vec::as_slice)
    }
}

impl fmt::Debug for Regions {
    /// How many regions there are and on how many chromosomes, and not
    /// each of them: a BED can hold a million.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Regions")
            .field("num_regions", &self.num_regions)
            .field("num_chroms", &self.of_each_chrom.len())
            .finish()
    }
}

/// Whether a line of a BED, without its end of line, is one that is
/// skipped: an empty one, a comment, and one whose first word, what comes
/// before the first space or tab, is `track` or `browser`. A chromosome
/// named `tracks1` is a region.
fn is_not_a_region(line: &[u8]) -> bool {
    let first_word = line
        .split(|byte| *byte == b' ' || *byte == b'\t')
        .next()
        .unwrap_or_default();
    line.is_empty() || line.starts_with(b"#") || WORDS_OF_THE_LINES_SKIPPED.contains(&first_word)
}

/// The chromosome and the region of one line of a BED that is not skipped,
/// or what is wrong with it.
fn region_of_the_line(line: &[u8]) -> std::result::Result<(&[u8], Region), BedLineProblem> {
    let mut columns = line.split(|byte| *byte == b'\t');
    let (Some(chrom), Some(start), Some(end)) = (columns.next(), columns.next(), columns.next())
    else {
        return Err(BedLineProblem::FewerThanThreeColumns {
            columns: line.split(|byte| *byte == b'\t').count(),
            with_spaces: line.contains(&b' '),
        });
    };
    if chrom.is_empty() {
        return Err(BedLineProblem::EmptyChromosome);
    }
    let start = whole_number(start).map_err(|problem| {
        let found = String::from_utf8_lossy(start).into_owned();
        match problem {
            NotANumberOf64Bits::NotDigits => BedLineProblem::StartNotAWholeNumber { found },
            NotANumberOf64Bits::AboveTheLargest => BedLineProblem::StartAboveTheLargest { found },
        }
    })?;
    let end = whole_number(end).map_err(|problem| {
        let found = String::from_utf8_lossy(end).into_owned();
        match problem {
            NotANumberOf64Bits::NotDigits => BedLineProblem::EndNotAWholeNumber { found },
            NotANumberOf64Bits::AboveTheLargest => BedLineProblem::EndAboveTheLargest { found },
        }
    })?;
    if start >= end {
        return Err(BedLineProblem::StartNotBelowEnd { start, end });
    }
    // The start is below the end, which is a `u64`, so it is below the
    // largest one and this does not saturate.
    let first = start.saturating_add(1);
    Ok((chrom, Region { first, last: end }))
}

/// Why a column of a BED is not a start or an end.
enum NotANumberOf64Bits {
    /// It is empty or holds a byte that is not a digit.
    NotDigits,
    /// It is digits alone, of a number above the largest `u64`.
    AboveTheLargest,
}

/// The whole number `digits` writes, or why it is not one that fits in a
/// `u64`.
fn whole_number(digits: &[u8]) -> std::result::Result<u64, NotANumberOf64Bits> {
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_digit) {
        return Err(NotANumberOf64Bits::NotDigits);
    }
    digits
        .iter()
        .try_fold(0_u64, |number, byte| {
            number
                .checked_mul(10)?
                .checked_add(u64::from(byte.saturating_sub(b'0')))
        })
        .ok_or(NotANumberOf64Bits::AboveTheLargest)
}

/// The regions sorted by their first position and joined where they
/// overlap or touch, in place: two regions touch when the first position of
/// one is the one after the last of the other.
///
/// The joined regions are written over the front of the vector, which is
/// then cut to them and given back the room it no longer needs, so no
/// second vector is built.
fn join(regions: &mut Vec<Region>) {
    regions.sort_unstable_by_key(|region| region.first);
    let mut num_joined: usize = 0;
    for next in 0..regions.len() {
        let Some(region) = regions.get(next).copied() else {
            break;
        };
        // The region the ones before `next` were joined into, which is at
        // `num_joined - 1` once there is one.
        let before = num_joined
            .checked_sub(1)
            .and_then(|last| regions.get_mut(last));
        match before {
            // The last position of a region is at most `u64::MAX`, and one
            // that ends there reaches every position after it, which the
            // saturated sum says too.
            Some(before) if region.first <= before.last.saturating_add(1) => {
                before.last = before.last.max(region.last);
            }
            Some(_) | None => {
                if let Some(slot) = regions.get_mut(num_joined) {
                    *slot = region;
                }
                // `num_joined` is at most `next`, which is below the length
                // of the vector, so this does not saturate.
                num_joined = num_joined.saturating_add(1);
            }
        }
    }
    regions.truncate(num_joined);
    regions.shrink_to_fit();
}

/// Whether `pos` is in one of `regions`, which are sorted and joined, by a
/// binary search: the first region that does not end before `pos` is the
/// only one that can hold it.
fn holds(regions: &[Region], pos: u64) -> bool {
    let at = regions.partition_point(|region| region.last < pos);
    regions.get(at).is_some_and(|region| region.first <= pos)
}

/// Which variants the filter by regions keeps: those inside the regions,
/// or, with `exclude`, those outside all of them.
///
/// It is what a reader is offered with [`BlockReader::skip_outside`], so
/// that a source that can pass over the variants the filter would take out
/// does not build them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionSelection {
    /// The regions, shared by the step and by every pass built from it.
    pub regions: Arc<Regions>,
    /// Whether the filter keeps the variants outside the regions and not
    /// those inside.
    pub exclude: bool,
}

impl RegionSelection {
    /// `"regions"`, or `"excluded_regions"` with `exclude`: the name the
    /// step and its counts have for a Python and a TypeScript user, and by
    /// which a second step of the same kind is refused.
    #[must_use]
    pub fn kind(&self) -> &'static str {
        if self.exclude {
            THE_KIND_OUTSIDE
        } else {
            THE_KIND_INSIDE
        }
    }

    /// Whether the variant at position `pos` of the chromosome named
    /// `chrom` is one the filter keeps.
    #[must_use]
    pub fn keeps(&self, chrom: &str, pos: u64) -> bool {
        self.regions_of(chrom.as_bytes()).keeps(pos)
    }

    /// Whether no position from `min_pos` to `max_pos` of `chrom`, both
    /// included, is one this selection keeps, which is what lets the vars
    /// file reader skip a batch. A `min_pos` above `max_pos` is no position
    /// at all, and none of them is kept.
    #[must_use]
    pub fn keeps_none_of(&self, chrom: &str, min_pos: u64, max_pos: u64) -> bool {
        self.regions_of(chrom.as_bytes())
            .keeps_none_of(min_pos, max_pos)
    }

    /// The selection on the chromosome named `chrom`, as the bytes a source
    /// holds the name in: a reader looks the name up once for a run of
    /// variants of one chromosome and asks the handle of each of them,
    /// where [`RegionSelection::keeps`] looks the name up at every call. A
    /// name the BED does not have is a chromosome of no region.
    #[must_use]
    pub fn regions_of(&self, chrom: &[u8]) -> SelectionOfAChrom<'_> {
        SelectionOfAChrom {
            regions: self.regions.of_the_chrom(chrom).unwrap_or(&[]),
            exclude: self.exclude,
        }
    }
}

/// A [`RegionSelection`] on one chromosome: its regions, sorted and joined,
/// and which side of them is kept.
#[derive(Debug, Clone, Copy)]
pub struct SelectionOfAChrom<'regions> {
    regions: &'regions [Region],
    exclude: bool,
}

impl SelectionOfAChrom<'_> {
    /// Whether the variant at position `pos` of this chromosome is one the
    /// filter keeps.
    #[must_use]
    pub fn keeps(&self, pos: u64) -> bool {
        holds(self.regions, pos) != self.exclude
    }

    /// Whether no position from `min_pos` to `max_pos` of this chromosome,
    /// both included, is one the filter keeps. A `min_pos` above `max_pos`
    /// is no position at all, and none of them is kept.
    #[must_use]
    pub fn keeps_none_of(&self, min_pos: u64, max_pos: u64) -> bool {
        if min_pos > max_pos {
            return true;
        }
        // The first region that does not end before `min_pos`: the regions
        // before it hold none of the positions, and it is the only one that
        // can hold `min_pos`.
        let first = self
            .regions
            .get(self.regions.partition_point(|region| region.last < min_pos))
            .copied();
        if self.exclude {
            // Every position is inside, and the regions are joined, so one
            // region holds them all.
            first.is_some_and(|region| region.first <= min_pos && max_pos <= region.last)
        } else {
            // No position is inside: the first region that could hold one
            // starts after the last of them.
            first.is_none_or(|region| region.first > max_pos)
        }
    }
}

/// The filter by regions of one pass: it keeps the variants of a block that
/// its selection keeps and counts how many it was given and how many it
/// kept.
///
/// It is an object of its own, apart from the reader that puts it over a
/// source, so that the rule and its counts are worked out on a block that
/// no reader gave, as [`VarFilter`](crate::filters::VarFilter) is.
#[derive(Debug)]
pub struct RegionFilter {
    selection: RegionSelection,
    stats: FilteringStats,
}

impl RegionFilter {
    /// The filter of `selection`, with both its counts at 0.
    #[must_use]
    pub fn new(selection: RegionSelection) -> RegionFilter {
        RegionFilter {
            selection,
            stats: FilteringStats::default(),
        }
    }

    /// Which regions it keeps the variants of, and on which side.
    #[must_use]
    pub fn selection(&self) -> &RegionSelection {
        &self.selection
    }

    /// The variants of the block that the selection keeps, kept in it in
    /// their order, and the others dropped from the genotypes and from
    /// every column, in place. `chroms` names the chromosome numbers of the
    /// block: it is the table of the reader the block came from.
    ///
    /// The variants of the block are added to the counts, and the ones
    /// kept to the ones kept. A block of no variants is left as it is.
    ///
    /// # Errors
    ///
    /// When the arrays of the block are not of the size the block states,
    /// which [`Block::check`] finds; when the block has variants and no
    /// chromosome or no position, which is the error of a field that is
    /// not in the block; and when a chromosome number of the block has no
    /// name in `chroms`, which a reader with a defect gives. After any of
    /// them the block is as it was and nothing was added to the counts.
    pub fn filter_block(&mut self, block: &mut Block, chroms: &ChromTable) -> Result<()> {
        block.check()?;
        let processed = block.num_vars;
        if processed == 0 {
            return Ok(());
        }
        let missing = Needs::CHROM_POS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        let (Some(chrom), Some(pos)) = (block.chrom.as_deref(), block.pos.as_deref()) else {
            // A block holds the chromosome and the position as one field
            // and only when both columns are there, so the check above is
            // what refuses a block without them and this is not reached.
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::CHROM_POS,
            });
        };
        let keep = keep_of_the_rows(&self.selection, chrom, pos, chroms)?;
        let kept = keep.iter().filter(|keep_it| **keep_it).count();
        block.retain_vars(&keep)?;
        // A `usize` is 64 bits on the targets popnei builds natively for
        // and 32 in wasm, so every one of them is a `u64` and neither
        // conversion takes the value it saturates at.
        self.stats.vars_processed = self
            .stats
            .vars_processed
            .saturating_add(u64::try_from(processed).unwrap_or(u64::MAX));
        self.stats.vars_kept = self
            .stats
            .vars_kept
            .saturating_add(u64::try_from(kept).unwrap_or(u64::MAX));
        Ok(())
    }

    /// How many variants it was given and how many it kept, over the
    /// blocks it was given. [`RegionsReader`] adds the variants its source
    /// passed over to the ones given when its counts are read.
    #[must_use]
    pub fn stats(&self) -> FilteringStats {
        self.stats
    }
}

/// The regions of the last chromosome a row was on, by its number, which a
/// thread keeps from one row to the next: the variants of a chromosome
/// come together in a sorted source, so the name is looked up once for
/// each run of them and not for each variant.
type TheLastChrom<'regions> = Option<(u32, SelectionOfAChrom<'regions>)>;

/// Whether the selection keeps the variant at `chrom_number` and `pos`,
/// with the regions of its chromosome read from `last` when the row before
/// was on the same one.
fn keeps_the_row<'regions>(
    selection: &'regions RegionSelection,
    chroms: &ChromTable,
    last: &mut TheLastChrom<'regions>,
    chrom_number: u32,
    pos: u64,
) -> Result<bool> {
    let of_the_chrom = match *last {
        Some((number, of_the_chrom)) if number == chrom_number => of_the_chrom,
        Some(_) | None => {
            let Some(name) = chroms.name(chrom_number) else {
                return Err(Error::RegionFilterChromNameMissing {
                    number: chrom_number,
                });
            };
            let of_the_chrom = selection.regions_of(name.as_bytes());
            *last = Some((chrom_number, of_the_chrom));
            of_the_chrom
        }
    };
    Ok(of_the_chrom.keeps(pos))
}

/// Which variants of a block the selection keeps, one value for each, in
/// the order of the block.
///
/// Natively the rows are read on the threads of rayon, as the threshold
/// filters read theirs: each row gives one value, collected in the order of
/// the block, so the values do not depend on how many threads there are.
///
/// # Errors
///
/// A chromosome number with no name in `chroms`: the error of the first
/// row that has one, wherever the threads found it, so that the same block
/// gives the same message every time.
#[cfg(not(target_family = "wasm"))]
fn keep_of_the_rows(
    selection: &RegionSelection,
    chrom: &[u32],
    pos: &[u64],
    chroms: &ChromTable,
) -> Result<Vec<bool>> {
    use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

    let keep: Result<Vec<bool>> = chrom
        .par_iter()
        .zip(pos.par_iter())
        .map_init(
            || None,
            |last, (chrom_number, position)| {
                keeps_the_row(selection, chroms, last, *chrom_number, *position)
            },
        )
        .collect();
    match keep {
        Ok(keep) => Ok(keep),
        // The rows are read again one after another to find the first one
        // that is an error, which costs a read of a block that is refused.
        Err(of_a_thread) => match keep_of_the_rows_one_by_one(selection, chrom, pos, chroms) {
            Err(of_the_first_row) => Err(of_the_first_row),
            Ok(_) => Err(of_a_thread),
        },
    }
}

/// The same values, with the rows read one after another, which is what
/// wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn keep_of_the_rows(
    selection: &RegionSelection,
    chrom: &[u32],
    pos: &[u64],
    chroms: &ChromTable,
) -> Result<Vec<bool>> {
    keep_of_the_rows_one_by_one(selection, chrom, pos, chroms)
}

/// The rows read one after another: what wasm does, and what the threads
/// fall back on to find the first row that is an error.
fn keep_of_the_rows_one_by_one(
    selection: &RegionSelection,
    chrom: &[u32],
    pos: &[u64],
    chroms: &ChromTable,
) -> Result<Vec<bool>> {
    let mut last = None;
    chrom
        .iter()
        .zip(pos)
        .map(|(chrom_number, position)| {
            keeps_the_row(selection, chroms, &mut last, *chrom_number, *position)
        })
        .collect()
}

/// A reader that gives the variants of its source that the filter by
/// regions keeps.
///
/// It offers the regions to its source when it is built, with
/// [`BlockReader::skip_outside`], and a source that takes them passes over
/// the variants the filter would take out, so that they are not built at
/// all. The filter still looks at every variant it is given, so the
/// variants it gives are the same whatever the source answered, and its
/// counts are too: the variants it was given are those that reached it and
/// those the source says it passed over.
///
/// It owns its source, since an offer the source took holds for as long as
/// the source lives. It keeps the contract of a reader of
/// `docs/specs/block.md`, as [`FilteredReader`](crate::filters::FilteredReader)
/// does: a block left with no variant is not given and the next one is
/// taken; after an error, of its source or of its filter, it gives `None`
/// at every call and does not call its source again; and a source that
/// gives a block of no variants has a defect and is the error of that.
pub struct RegionsReader<R: BlockReader> {
    reader: R,
    filter: RegionFilter,
    /// Whether the source took the regions, and so passes over the
    /// variants outside them and counts those in its `num_skipped`.
    source_skips: bool,
    /// Whether the source has no more blocks or one of the two, the source
    /// or the filter, gave an error. After any of them there is no block.
    finished: bool,
}

impl<R: BlockReader> RegionsReader<R> {
    /// The reader that gives the variants of `reader` that `filter` keeps.
    ///
    /// It offers the selection of `filter` to `reader`, which takes it or
    /// not; the filter reads every variant it is given either way. It asks
    /// its source for the chromosome and the position besides what its
    /// consumer asks for, when the consumer calls
    /// [`BlockReader::set_needs`].
    ///
    /// # Errors
    ///
    /// When `reader` holds a filter by regions of the kind of `filter`
    /// already, which its [`BlockReader::filtering_stats`] says: two sets of
    /// regions on one side are one set, which one BED says. No offer is
    /// made then.
    pub fn new(mut reader: R, filter: RegionFilter) -> Result<RegionsReader<R>> {
        let kind = filter.selection().kind();
        if reader
            .filtering_stats()
            .iter()
            .any(|(of_the_chain, _)| *of_the_chain == kind)
        {
            return Err(Error::RegionFilterOfAKindThatIsSet { kind });
        }
        let source_skips = reader.skip_outside(filter.selection().clone());
        Ok(RegionsReader {
            reader,
            filter,
            source_skips,
            finished: false,
        })
    }
}

impl<R: BlockReader> BlockReader for RegionsReader<R> {
    /// The next block of the source with the variants the selection keeps
    /// kept in it, and the blocks its filter emptied passed over.
    ///
    /// # Errors
    ///
    /// When the source fails; when a block of the source holds no variant,
    /// which no reader of popnei gives; and what
    /// [`RegionFilter::filter_block`] refuses: a block whose arrays are not
    /// of its size, a block that has variants and no chromosome or no
    /// position, and a chromosome number the table of the source has no
    /// name for. After any of them there is no block and the source is not
    /// called again.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        loop {
            let mut block = match self.reader.next_block() {
                Ok(Some(block)) => block,
                Ok(None) => {
                    self.finished = true;
                    return Ok(None);
                }
                Err(error) => {
                    self.finished = true;
                    return Err(error);
                }
            };
            // A source that gives a block of no variants has a defect, and
            // it is not asked again: over a source that always gives one, a
            // reader that asked again would never come back.
            if block.num_vars == 0 {
                self.finished = true;
                return Err(Error::ReaderGaveABlockOfNoVariants);
            }
            if let Err(error) = self.filter.filter_block(&mut block, self.reader.chroms()) {
                self.finished = true;
                return Err(error);
            }
            if block.num_vars > 0 {
                return Ok(Some(block));
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

    /// The fields of the consumer, the chromosome and the position, which
    /// the filter reads for every variant: so the blocks this reader gives
    /// hold both also when the consumer asked for neither.
    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs.union(Needs::CHROM_POS));
    }

    /// The counts of this filter, with the variants its source passed over
    /// among the ones given when the source took the regions, and after
    /// them those of the filters between the source and its own source.
    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        let mut own = self.filter.stats();
        if self.source_skips {
            own.vars_processed = own.vars_processed.saturating_add(self.reader.num_skipped());
        }
        let mut stats = vec![(self.filter.selection().kind(), own)];
        stats.extend(self.reader.filtering_stats());
        stats
    }

    fn header(&self) -> &SourceHeader {
        self.reader.header()
    }

    /// False, and the source is not offered the regions: the variants it
    /// passed over would never reach this filter, and its counts would be
    /// of fewer variants than the source gave.
    fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
        false
    }

    /// 0: the variants its own source passed over are counted in its own
    /// counts, and a filter by regions over this one adds nothing to them.
    fn num_skipped(&self) -> u64 {
        0
    }
}

impl<R: BlockReader> fmt::Debug for RegionsReader<R> {
    /// What it filters and where it has got to. The source is left out, so
    /// that a `RegionsReader` over a reader that has no `Debug` has one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("RegionsReader")
            .field("filter", &self.filter)
            .field("source_skips", &self.source_skips)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests;
