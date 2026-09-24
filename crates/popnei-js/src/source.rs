//! What the two sources of variants share: one pass over a source, what that
//! pass tells the page while it reads, the columns of the blocks it gives
//! JavaScript, the bytes of the vars file written from it, and how many
//! passes each consumer of the package makes.
//!
//! A source is a file with what is needed to read it, a VCF with its
//! options in `vcf.rs` and a vars file in `vars.rs`, and the file is in one
//! of the two places of [`TheFileOfASource`]: a copy of the whole of it in
//! the memory of wasm, the bytes of a `Uint8Array` the application gave, or
//! a file the user picked in the page, which is read one range at a time and
//! is never in that memory whole. Each source is a
//! class of wasm-bindgen that the `Variants` of the TypeScript package
//! holds, each reads what says what the file holds when it is built, so
//! bytes that are not of its format fail at `openVcf` or at `openVars`, and
//! each reads the bytes again at every pass, which is what lets a user give
//! the same `Variants` to one calculation after another. What they have in
//! common is [`OpenSource`], the reader of one pass.
//!
//! [`PassOverTheBytes`] is where the bytes of a pass come from, and it is the
//! one part of a pass that comes back out to JavaScript: it counts what the
//! pass has read and tells the page at its first read and once per range of
//! bytes after that, so that a page can draw a bar over a run. The end of a
//! run tells the page once more for each of its passes, which is what says
//! that a pass is over, because no read does. A run is one call of one
//! consumer with the passes it makes, [`Run`] in [`RUNS`], and what a source
//! keeps in
//! JavaScript, the function it tells, is [`InJavaScript`] in
//! [`IN_JAVASCRIPT`]: no handle of JavaScript is `Send`, and a reader of the
//! core has to be, so what a pass holds of those two tables is the number of
//! an entry. A function of the page that throws ends the pass where it was
//! reading, and the value it threw is kept in the run and given back to the
//! application in place of the error the failed read became, so a cancel is
//! the application's own value and popnei's errors stay the ones it made
//! itself.
//!
//! [`Blocks`] is that pass, whichever source it came from: it owns the chain
//! of readers of the pass, the source with a filter over it for each step of
//! the `Variants` and a `Reblock` at its end, which gives the blocks the
//! size that was asked for. It counts the variants of the blocks it gives
//! and reads the counts of the filters from that chain, which is the
//! [`PassCounts`] that the package turns into the `passStats` of
//! `docs/specs/variant.md`. [`BlockColumns`] is one block on its way out,
//! and [`VarsWritten`] the bytes of a vars file with the counts of the pass
//! that wrote it.
//!
//! They live in the memory of wasm, which the garbage collector of
//! JavaScript does not see, so the TypeScript package frees each of them:
//! the block as soon as its columns are read, the counts as soon as their
//! numbers are read, the pass when the iteration ends or is abandoned, and
//! the source when the user calls `free()`.
//!
//! The columns leave as JavaScript arrays, which wasm-bindgen copies out of
//! the memory of wasm: a typed array that is a view into that memory stops
//! being valid when it grows, as section 11 of `docs/architecture.md` says.
//! Each column is moved out of the block as it is read, so that the copy
//! that crosses is the only one.

use std::cell::RefCell;
use std::io::{BufRead, Cursor, ErrorKind, Read, Seek, SeekFrom, Write};
use std::sync::Arc;

use js_sys::{Array, Function, Uint8Array};
use wasm_bindgen::prelude::wasm_bindgen;
use wasm_bindgen::{JsCast, JsValue};
use web_sys::{Blob, FileReaderSync};

use popnei::block::{AllelesColumn, Block, BlockReader, Reblock, needs_of_the_fields};
use popnei::filters::FilteringStats;
use popnei::variant::Needs;

use crate::errors::JsPopneiError;
use crate::steps::{Steps, chain_of};

/// The largest position a block hands to JavaScript, 2^53, which
/// `docs/specs/block.md` sets.
///
/// The positions cross as float64, which holds every whole number up to
/// this one and not the ones above it: 2^53 + 1 would arrive as 2^53. It is
/// also the longest stretch of a chromosome the resampling groups of the
/// distances between populations are cut into, for the same reason.
///
/// It is one above the largest window a user may write for the filter by
/// linkage disequilibrium, the 2^53 - 1 of `LARGEST_WINDOW` of `steps.rs`,
/// and the two are not the same kind of number. A position is read from a
/// file, and 2^53 itself is held exactly, so it is handed out. A window is
/// written by a user and read back to them, which is what
/// `Number.isSafeInteger` stands for: above 2^53 - 1 the numbers a user
/// can write no longer run one by one.
pub(crate) const LARGEST_POSITION: u64 = 9_007_199_254_740_992;

/// One consumer of the package, the function that runs a `Variants`, with
/// the argument of the one whose number of passes depends on it.
///
/// One call of a consumer is a run, and a run is one or two passes over the
/// source, each a reading of it from its start, which is how
/// `docs/glossary.md` has the three words. Which consumer it is says how
/// many passes, and nothing else does, so this is what a consumer names
/// itself with when it asks for [`Consumer::num_passes`].
pub(crate) enum Consumer {
    /// `calcPerVarDistribs`.
    PerVarDistribs,
    /// `calcPerIndividualStats`.
    PerIndividualStats,
    /// `calcPairwiseKosmanDists`.
    KosmanDists,
    /// `calcPopDists`.
    PopDists,
    /// `calcRogersHuffR2Matrix`.
    R2Matrix,
    /// `calcKinship`.
    Kinship,
    /// `doPcaFromVariants`, the one consumer that reads the source twice.
    PcaOfVariants {
        /// How many components the weight of each variant is asked for,
        /// the `numPrinComps` of the call, where 0 asks for no weight.
        num_prin_comps: usize,
    },
    /// `calcGwas`.
    Gwas,
    /// `writeVars`.
    WriteVars,
    /// The iteration of `iterBlocks`.
    IterBlocks,
}

impl Consumer {
    /// How many passes over the source this consumer makes.
    ///
    /// Every consumer reads the source once, except the principal
    /// components of the variants asked for weights: a weight needs the
    /// eigenvectors, which are known when the first pass ends, so the
    /// variants are read a second time. `numPrinComps` 0 asks for no
    /// weight and reads the source once.
    ///
    /// It is what [`num_passes_of`] gives a page before a run starts, and
    /// it is also what a consumer opens its run with, so that every call
    /// that tells the page how far a pass has got carries this same number
    /// and the two cannot disagree.
    pub(crate) fn num_passes(&self) -> u32 {
        match *self {
            Consumer::PcaOfVariants { num_prin_comps } => {
                if num_prin_comps > 0 {
                    2
                } else {
                    1
                }
            }
            Consumer::PerVarDistribs
            | Consumer::PerIndividualStats
            | Consumer::KosmanDists
            | Consumer::PopDists
            | Consumer::R2Matrix
            | Consumer::Kinship
            | Consumer::Gwas
            | Consumer::WriteVars
            | Consumer::IterBlocks => 1,
        }
    }

    /// The consumer a user of the package named, with the `num_prin_comps`
    /// of the call, which every consumer but the principal components of
    /// the variants ignores.
    ///
    /// # Errors
    ///
    /// When `name` is of no consumer of the package.
    fn of_the_name(name: &str, num_prin_comps: usize) -> Result<Consumer, JsPopneiError> {
        match name {
            "calcPerVarDistribs" => Ok(Consumer::PerVarDistribs),
            "calcPerIndividualStats" => Ok(Consumer::PerIndividualStats),
            "calcPairwiseKosmanDists" => Ok(Consumer::KosmanDists),
            "calcPopDists" => Ok(Consumer::PopDists),
            "calcRogersHuffR2Matrix" => Ok(Consumer::R2Matrix),
            "calcKinship" => Ok(Consumer::Kinship),
            "doPcaFromVariants" => Ok(Consumer::PcaOfVariants { num_prin_comps }),
            "calcGwas" => Ok(Consumer::Gwas),
            "writeVars" => Ok(Consumer::WriteVars),
            "iterBlocks" => Ok(Consumer::IterBlocks),
            _ => Err(JsPopneiError::Refused(format!(
                "`{name}` is not a consumer of popnei, which are the functions \
                 that read the variants of a `Variants`: {names}",
                names = THE_CONSUMERS.join(", ")
            ))),
        }
    }
}

/// The name of each consumer as a user of the package writes it, for the
/// message of a name that is of none of them.
const THE_CONSUMERS: [&str; 10] = [
    "calcPerVarDistribs",
    "calcPerIndividualStats",
    "calcPairwiseKosmanDists",
    "calcPopDists",
    "calcRogersHuffR2Matrix",
    "calcKinship",
    "doPcaFromVariants",
    "calcGwas",
    "writeVars",
    "iterBlocks",
];

/// How many passes over the source the consumer called `consumer` makes,
/// with the `num_prin_comps` of the call, which every consumer but
/// `doPcaFromVariants` ignores.
///
/// A page that draws one bar for a whole run asks this before the run
/// starts, and every call that tells the page how far a pass has got
/// carries this same number: the consumer asks for it here and opens its
/// run with it, so the two cannot disagree.
///
/// # Errors
///
/// When `consumer` is of no consumer of the package.
#[wasm_bindgen]
pub fn num_passes_of(consumer: &str, num_prin_comps: usize) -> Result<u32, JsPopneiError> {
    Ok(Consumer::of_the_name(consumer, num_prin_comps)?.num_passes())
}

/// A file of variants that was opened, which every pass reads again.
pub(crate) trait OpenSource {
    /// How many alleles the genotype of one individual holds in the blocks
    /// of every pass over the source, which a calculation that counts
    /// genotypes out of called alleles is built with.
    fn ploidy(&self) -> usize;

    /// The run of `consumer` over this source, which every reader the
    /// consumer opens belongs to and which is taken out of [`RUNS`] when it
    /// is dropped.
    ///
    /// A consumer opens it before its first reader and holds it until it has
    /// its result, and `iterBlocks` holds it in the [`Blocks`] of its
    /// iteration: a pass takes its number from the run when it first reads,
    /// and the run is what says how many passes the page is told the run
    /// makes.
    fn starts_a_run(&self, consumer: &Consumer) -> RunOfAConsumer;

    /// The reader of one pass of `run` over the source, which reads the bytes
    /// from their start.
    ///
    /// `num_vars_per_block` is the size the caller will ask the blocks for,
    /// which a source whose reader can give them at that size is built
    /// with, so that the `Reblock` over it has nothing to cut or to join. A
    /// source that gives its blocks as they are in the file, the vars file
    /// whose batches were written at one size, leaves it to that `Reblock`.
    ///
    /// # Errors
    ///
    /// When what says what the file holds, the header of a VCF or the
    /// schema of a vars file, cannot be read.
    fn reader(
        &self,
        run: &RunOfAConsumer,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error>;
}

/// The bytes of a file, which every pass over it shares.
///
/// `Cursor` reads anything that gives a `&[u8]`, and an `Arc<Vec<u8>>` gives
/// a `&Vec<u8>`, so this is the one line that turns the second into the
/// first. What it holds is the `Vec` wasm-bindgen filled with the bytes of
/// the `Uint8Array`, and nothing copies it again: a copy of a file of 80 MB
/// stays in the memory of wasm for as long as the tab lives, because that
/// memory grows and never shrinks.
pub(crate) struct SharedBytes(Arc<Vec<u8>>);

impl AsRef<[u8]> for SharedBytes {
    fn as_ref(&self) -> &[u8] {
        self.0.as_slice()
    }
}

/// The bytes of a source read from their first one, which every pass over
/// that source shares.
fn cursor_of(bytes: &Arc<Vec<u8>>) -> Cursor<SharedBytes> {
    Cursor::new(SharedBytes(Arc::clone(bytes)))
}

/// How many bytes a pass over a file of the page reads at a time, and how
/// many bytes any pass reads between two calls that tell the page how far it
/// has got, 4 MiB.
///
/// The two are one number. A pass over a file of the page is told of its
/// progress once per range it reads, which is the only moment where it has
/// read a known number of bytes more, and a pass over bytes that are already
/// in the memory of wasm tells the page as often as that one does, so that a
/// bar moves the same way over both.
///
/// Nothing has been measured at this number: "Speed" of
/// `docs/specs/js_sources.md` leaves it at 4 MiB until work package 4 of
/// `docs/plans/js-sources.md` times one pass over a VCF of a few hundred MB
/// in Chromium at 256 KiB, 1 MiB, 4 MiB and 16 MiB, and sets it from what it
/// measures.
const NUM_BYTES_PER_RANGE: u64 = 4 * 1024 * 1024;

/// The number of no entry of [`RUNS`] or of [`IN_JAVASCRIPT`], which a run
/// that could not be put in the first carries: its passes are told to nobody.
///
/// No tab reaches it. An entry of `RUNS` is there while its run is, and a run
/// holds the readers of its passes in the memory of wasm, which is 4 GB.
const NO_ENTRY: u32 = u32::MAX;

/// One pass over the bytes of a source, which the readers of the core take as
/// they take a cursor over an array of bytes.
///
/// It counts what the pass has read and tells the page once per range of
/// bytes, and it is `Send`, because what it holds of JavaScript is the number
/// of an entry of [`RUNS`] and not a handle. `BlockReader`, which everything
/// that gives blocks implements, asks for `Send`, as section 1 of
/// `docs/architecture.md` says.
///
/// What `Send` means here is that the pass moves and is then told nothing:
/// [`RUNS`] and [`IN_JAVASCRIPT`] are of the thread that made the source, so
/// a pass that is read from another thread finds neither its run nor the
/// function of the page, takes no number, tells nobody and gives the same
/// bytes as ever. A tab has one thread of wasm, and moving a reader is what
/// the read ahead thread of section 3 of `docs/architecture.md` would do
/// natively.
pub(crate) struct PassOverTheBytes {
    bytes: TheBytes,
    /// Which run of [`RUNS`] this pass belongs to, which is what says which
    /// source it reads, which pass of the run it is and how many there are.
    run: u32,
    /// Which pass of the run this is, 1 for the first, and 0 while
    /// `took_its_number` is false and for a pass whose run is not in the
    /// table of this thread.
    pass: u32,
    /// Whether the first read has asked the run for the number of this pass,
    /// which happens once and not at every read of a pass that got 0: the
    /// number comes from a table of the thread that made the source, and a
    /// pass that was moved to another thread finds none there.
    ///
    /// The number is taken at the first read and not when the pass is built,
    /// because the principal components of the variants build both of their
    /// readers before either of them gives a block.
    took_its_number: bool,
    /// How many bytes this pass has read, which is never more than
    /// `num_bytes`: a pass over a vars file reads its footer and its batches
    /// and not the whole of it.
    bytes_read: u64,
    /// How many it had read when the page was last told, which is what the
    /// size of a range is compared against to decide whether to tell it
    /// again.
    told_at: u64,
    num_bytes: u64,
    /// Why this pass gives an error instead of bytes, and nothing while it
    /// reads.
    ended: Option<WhyThePassEnded>,
}

/// Why a pass gives an error instead of bytes before the source is over.
///
/// Every read after the one it happened in fails with the same error, and
/// the page is told no more of that pass, the call that would say it ended
/// among them, as `docs/specs/js_sources.md` says. No reader of the core
/// reads again after an error of either kind, so what this keeps is the
/// promise made to the application and not the behaviour of any of them.
enum WhyThePassEnded {
    /// The function of the page threw, which is how an application stops a
    /// run.
    TheApplicationStopped,
    /// A read gave more bytes than a count of 64 bits holds, which is a
    /// defect of this crate: no target popnei builds for has a `usize` wider
    /// than a `u64`.
    TheCountDidNotFit,
}

/// Where the bytes of a pass come from.
enum TheBytes {
    /// A copy of the whole file in the memory of wasm, which every pass over
    /// that source shares.
    InMemory(Cursor<SharedBytes>),
    /// A file of the page, read one range at a time.
    OfAFile(RangesOfAFile),
}

/// One pass over a file of the page: the range of bytes it holds, where that
/// range starts in the file, and where the pass is.
///
/// A read gives what is left of the range, and a read that starts where the
/// range ends asks the browser for the next one, [`NUM_BYTES_PER_RANGE`]
/// bytes or what is left of the file, whichever is fewer. So what a pass over
/// a file holds in the memory of wasm is one range and what the reader over
/// it builds, and the file is never there whole.
///
/// A seek moves where the pass is and reads nothing: the range is read again
/// only when the pass reads outside the one it holds, which is what lets the
/// reader of a vars file jump to its footer and back inside one range with no
/// second call into JavaScript.
///
/// The `Blob` and the `FileReaderSync` are not here. A reader of the core has
/// to be `Send` and no handle of JavaScript is, so they live in the entry of
/// [`IN_JAVASCRIPT`] that `source` numbers, which no thread leaves.
struct RangesOfAFile {
    /// Which entry of [`IN_JAVASCRIPT`] holds the file this pass reads and
    /// the reader of its ranges.
    source: u32,
    /// The bytes of the range the pass holds, which are none until its first
    /// read.
    range: Vec<u8>,
    /// Where that range starts in the file.
    range_at: u64,
    /// Where the pass reads next, which a seek moves.
    pos: u64,
    /// How many bytes the file holds, which `Blob.size` gave when the source
    /// was opened.
    ///
    /// It is the `num_bytes` of the [`PassOverTheBytes`] this is the bytes
    /// of, put here as well because the end of the file is what a read and a
    /// seek from the end are against. Neither is written after the pass is
    /// built, so the two cannot come apart.
    num_bytes: u64,
}

impl RangesOfAFile {
    /// The bytes of the range the pass holds from where the pass is, which
    /// are none when the pass is outside that range.
    fn what_it_holds(&self) -> &[u8] {
        let Some(from) = self.pos.checked_sub(self.range_at) else {
            return &[];
        };
        let Ok(from) = usize::try_from(from) else {
            return &[];
        };
        self.range.get(from..).unwrap_or(&[])
    }

    /// The bytes the file has ready, which are what is left of the range the
    /// pass holds, and the range that starts where the pass is when it holds
    /// none of it.
    ///
    /// One range is read from the browser and no more: what the pass holds is
    /// looked at with the arithmetic of `what_it_holds`, so the end of a
    /// range costs one call into JavaScript and not two.
    ///
    /// # Errors
    ///
    /// Those of [`RangesOfAFile::reads_the_range`]: a `Blob` that gave no
    /// range or no bytes, and a range that came back short.
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        if self.pos >= self.num_bytes {
            return Ok(&[]);
        }
        if self.what_it_holds().is_empty() {
            self.reads_the_range()?;
        }
        Ok(self.what_it_holds())
    }

    /// Reads the range of the file that starts where the pass is and leaves
    /// it as the range the pass holds.
    ///
    /// The `Blob` and the reader of its ranges are cloned out of
    /// [`IN_JAVASCRIPT`], which is two handles of JavaScript copied, so that
    /// no table of this crate is borrowed while the browser reads.
    ///
    /// # Errors
    ///
    /// When the source is not in the table of this thread, which is a pass
    /// that was moved to another thread or a source that was freed under it;
    /// when `Blob.slice` or `FileReaderSync` throws, a file that was moved or
    /// changed on disk among the causes; and when the range comes back
    /// shorter than the one that was asked for, which is not the end of the
    /// file.
    fn reads_the_range(&mut self) -> std::io::Result<()> {
        let num_bytes = NUM_BYTES_PER_RANGE.min(self.num_bytes.saturating_sub(self.pos));
        let at = self.pos;
        let file =
            IN_JAVASCRIPT.with_borrow(|sources| entry_of(sources, self.source)?.file.clone());
        let Some((blob, reader)) = file else {
            return Err(std::io::Error::other(format!(
                "the file of this pass is not among the sources popnei has open, so the \
                 {num_bytes} bytes from {at} of it cannot be read: a file of the page is \
                 read from the thread its source was opened in and from no other"
            )));
        };
        let Some(end) = at.checked_add(num_bytes) else {
            return Err(std::io::Error::other(format!(
                "the {num_bytes} bytes from {at} of this file end beyond what a count of \
                 64 bits holds, which is a defect of popnei; please report it"
            )));
        };
        // The two ends of the range cross as float64, which holds every whole
        // number up to 2^53. The size of the file was checked against that
        // number when the source was opened, so every position inside it is
        // the number popnei asked for.
        let range = blob
            .slice_with_f64_and_f64(at as f64, end as f64)
            .map_err(|thrown| self.the_browser_refused(at, num_bytes, &thrown))?;
        let buffer = reader
            .read_as_array_buffer(&range)
            .map_err(|thrown| self.the_browser_refused(at, num_bytes, &thrown))?;
        let bytes = Uint8Array::new(&buffer);
        let num_given = u64::from(bytes.length());
        if num_given != num_bytes {
            return Err(std::io::Error::other(format!(
                "popnei asked this file for the {num_bytes} bytes from {at} and the \
                 browser gave {num_given} of them, in a file of {num_bytes_of_the_file} \
                 bytes: a range that comes back short inside a file of that size is a \
                 file that changed after the page got its handle, and popnei ends the \
                 pass here instead of reading it as the end of the file. Pick the file \
                 again.",
                num_bytes_of_the_file = self.num_bytes
            )));
        }
        self.range = bytes.to_vec();
        self.range_at = at;
        Ok(())
    }

    /// What a range the browser refused fails with: what it threw, with the
    /// range that was asked for and the size of the file.
    fn the_browser_refused(&self, at: u64, num_bytes: u64, thrown: &JsValue) -> std::io::Error {
        std::io::Error::other(format!(
            "the browser did not give popnei the {num_bytes} bytes from {at} of this \
             file, which holds {num_bytes_of_the_file} bytes, and said: {said}",
            num_bytes_of_the_file = self.num_bytes,
            said = what_javascript_said(thrown)
        ))
    }

    /// The bytes of the file from where the pass is into `buf`, at most what
    /// is left of the range it holds.
    ///
    /// # Errors
    ///
    /// Those of [`RangesOfAFile::fill_buf`].
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let num_read = {
            let ready = self.fill_buf()?;
            let num_read = ready.len().min(buf.len());
            let (Some(into), Some(from)) = (buf.get_mut(..num_read), ready.get(..num_read)) else {
                return Err(std::io::Error::other(
                    "a read of a file of the page could not take the bytes of the range \
                     it holds, which is a defect of popnei; please report it",
                ));
            };
            into.copy_from_slice(from);
            num_read
        };
        self.consume(num_read);
        Ok(num_read)
    }

    /// The pass read `num_bytes` of the range it holds.
    ///
    /// Where the pass is moves by that many and is not held at the end of
    /// the file, which is what a `Cursor` over bytes in memory does with the
    /// same call: a pass that was moved beyond the end of the file keeps the
    /// position it was moved to, and its reads give no byte.
    fn consume(&mut self, num_bytes: usize) {
        let num_bytes = u64::try_from(num_bytes).unwrap_or(u64::MAX);
        self.pos = self.pos.saturating_add(num_bytes);
    }

    /// Where the pass reads next, which reads no byte: the range is read when
    /// the pass reads outside the one it holds.
    ///
    /// # Errors
    ///
    /// When the position asked for is before the first byte of the file or
    /// beyond what a count of 64 bits holds, which is what a `Cursor` over
    /// bytes in memory refuses there too.
    fn seek(&mut self, to: SeekFrom) -> std::io::Result<u64> {
        let pos = match to {
            SeekFrom::Start(at) => Some(at),
            SeekFrom::End(at) => self.num_bytes.checked_add_signed(at),
            SeekFrom::Current(at) => self.pos.checked_add_signed(at),
        };
        let Some(pos) = pos else {
            return Err(std::io::Error::new(
                ErrorKind::InvalidInput,
                "the pass was moved to a position that is before the first byte of the \
                 file or beyond what a count of 64 bits holds",
            ));
        };
        self.pos = pos;
        Ok(pos)
    }
}

/// What JavaScript threw, for the message of a range of a file that could not
/// be read: the text of the value when it is one, the `message` of an
/// `Error` when it is one, and what `Debug` writes of it for anything else.
///
/// A browser that refuses a range throws a `DOMException`, whose prototype
/// chain holds `Error`, so what a `NotReadableError` gives here is the
/// sentence it carries and not the `JsValue(...)` of the `Debug`.
///
/// It is the browser's sentence and not popnei's, and the message that
/// carries it says so.
fn what_javascript_said(thrown: &JsValue) -> String {
    if let Some(text) = thrown.as_string() {
        return text;
    }
    if let Some(error) = thrown.dyn_ref::<js_sys::Error>() {
        return String::from(error.message());
    }
    format!("{thrown:?}")
}

/// Where the file of a source is, which every pass over it reads again from
/// its first byte.
pub(crate) enum TheFileOfASource {
    /// A copy of the whole file in the memory of wasm, the `Vec`
    /// wasm-bindgen filled with the bytes of a `Uint8Array`, which every pass
    /// over the source shares.
    InMemory(Arc<Vec<u8>>),
    /// A file the user picked in the page, read one range at a time through
    /// the `Blob` and the `FileReaderSync` of the entry of [`IN_JAVASCRIPT`]
    /// that the source holds the number of, with the bytes `Blob.size` said
    /// it holds.
    OfThePage {
        /// How many bytes the file holds.
        num_bytes: u64,
    },
}

impl TheFileOfASource {
    /// The reader of what says what the file holds, the header of a VCF or
    /// the schema and the footer of a vars file, which `openVcf` and
    /// `openVars` read before they return.
    ///
    /// `source` is the number of the entry of [`IN_JAVASCRIPT`] the source
    /// keeps its file in. The pass belongs to no run, because there is no
    /// `Variants` yet to start one over, so it takes no number and is told to
    /// nobody; every reading that a consumer makes goes through
    /// [`TheFileOfASource::a_pass_of`] instead.
    ///
    /// # Errors
    ///
    /// Those of [`TheFileOfASource::a_pass_of`].
    pub(crate) fn the_opening_pass(&self, source: u32) -> Result<PassOverTheBytes, popnei::Error> {
        self.a_pass(source, NO_ENTRY)
    }

    /// One pass of `run` over the file, from its first byte.
    ///
    /// # Errors
    ///
    /// When the bytes in memory are more than a `u64` counts, which no target
    /// popnei builds for reaches: a `usize` is 32 bits in wasm and 64
    /// natively. The size of the file is what every call that tells the page
    /// carries, so a conversion that could not be made is an error of this
    /// crate and not a `numBytes` of 18446744073709551615.
    pub(crate) fn a_pass_of(
        &self,
        source: u32,
        run: &RunOfAConsumer,
    ) -> Result<PassOverTheBytes, popnei::Error> {
        self.a_pass(source, run.0)
    }

    /// One pass of the run numbered `run` over the file, from its first byte.
    ///
    /// # Errors
    ///
    /// Those of [`TheFileOfASource::a_pass_of`].
    fn a_pass(&self, source: u32, run: u32) -> Result<PassOverTheBytes, popnei::Error> {
        let (bytes, num_bytes) = match *self {
            TheFileOfASource::InMemory(ref bytes) => {
                let num_bytes = u64::try_from(bytes.len()).map_err(|_| {
                    popnei::Error::Io(std::io::Error::other(format!(
                        "the source holds {num_bytes} bytes, more than the count of a \
                         pass over it holds",
                        num_bytes = bytes.len()
                    )))
                })?;
                (TheBytes::InMemory(cursor_of(bytes)), num_bytes)
            }
            TheFileOfASource::OfThePage { num_bytes } => (
                TheBytes::OfAFile(RangesOfAFile {
                    source,
                    range: Vec::new(),
                    range_at: 0,
                    pos: 0,
                    num_bytes,
                }),
                num_bytes,
            ),
        };
        Ok(PassOverTheBytes {
            bytes,
            run,
            pass: 0,
            took_its_number: false,
            bytes_read: 0,
            told_at: 0,
            num_bytes,
            ended: None,
        })
    }
}

impl PassOverTheBytes {
    /// What a read of the pass does before it takes bytes: the first read
    /// takes the number of the pass from its run and tells the page that it
    /// has read nothing, and a later read tells it again when the pass has
    /// read a range of bytes since the last call.
    ///
    /// # Errors
    ///
    /// When the function of the page throws, and every read after the one it
    /// threw in.
    fn before_a_read(&mut self) -> std::io::Result<()> {
        match self.ended {
            Some(WhyThePassEnded::TheApplicationStopped) => {
                return Err(the_pass_was_stopped());
            }
            Some(WhyThePassEnded::TheCountDidNotFit) => {
                return Err(the_count_did_not_fit());
            }
            None => {}
        }
        if !self.took_its_number {
            self.took_its_number = true;
            self.pass = the_pass_that_starts(self.run);
            return self.tell();
        }
        if self.bytes_read.saturating_sub(self.told_at) >= NUM_BYTES_PER_RANGE {
            return self.tell();
        }
        Ok(())
    }

    /// The `num_read` bytes a read of the pass gave, counted against the size
    /// of the file.
    ///
    /// # Errors
    ///
    /// When `num_read` is more than a `u64` counts, which no target popnei
    /// builds for reaches. The pass ends there: a count that was not made is
    /// a bar that stands still, and every read after it fails with the same
    /// error.
    fn has_read(&mut self, num_read: usize) -> std::io::Result<()> {
        let Ok(num_read) = u64::try_from(num_read) else {
            self.ended = Some(WhyThePassEnded::TheCountDidNotFit);
            return Err(the_count_did_not_fit());
        };
        self.bytes_read = self.bytes_read.saturating_add(num_read).min(self.num_bytes);
        Ok(())
    }

    /// Tells the page how far this pass has got, and does nothing when the
    /// source of the run was given no function.
    ///
    /// The function is called with no table of this crate borrowed, so an
    /// application that calls popnei from inside it does not trap. What it
    /// throws is kept in the run, which is where the consumer of that run
    /// reads it to throw it in place of the error this read fails with.
    ///
    /// # Errors
    ///
    /// When the function throws, which ends the pass where it was reading.
    fn tell(&mut self) -> std::io::Result<()> {
        self.told_at = self.bytes_read;
        let Some((told, num_passes)) = what_tells_the_page(self.run) else {
            return Ok(());
        };
        // The four numbers of the `Progress` of `docs/specs/js_sources.md`,
        // which the package puts in the object its user reads: how many bytes
        // of the file this pass has read, how many the file holds, which pass
        // of the run is reading and how many passes the run makes. A count of
        // bytes of a file that is in the memory of a tab is below 2^32 and is
        // exact as a float64.
        let progress = Array::of4(
            &JsValue::from_f64(self.bytes_read as f64),
            &JsValue::from_f64(self.num_bytes as f64),
            &JsValue::from_f64(f64::from(self.pass)),
            &JsValue::from_f64(f64::from(num_passes)),
        );
        if let Err(thrown) = told.apply(&JsValue::NULL, &progress) {
            // The value is put in the run before the read fails, so that the
            // consumer finds it there whatever the readers of the core make
            // of the failed read on the way out.
            the_run_was_stopped(self.run, thrown);
            self.ended = Some(WhyThePassEnded::TheApplicationStopped);
            return Err(the_pass_was_stopped());
        }
        Ok(())
    }
}

impl Drop for PassOverTheBytes {
    /// The pass is over, which its run keeps so that its end can tell the
    /// page how far this pass got.
    ///
    /// No read says that a pass is over: a pass over a vars file stops after
    /// its last batch, with up to a range of bytes read since the last call,
    /// and a run that fails stops wherever it failed. So the reader of the
    /// pass being dropped is what says it, and the run makes the call when
    /// it is dropped in its turn, after every reader of it.
    ///
    /// A pass that never read is not among them, which is a pass whose
    /// reader was never built, since a reader reads when it is built; and
    /// neither is a pass the application stopped, which is not told how far
    /// it had got.
    fn drop(&mut self) {
        if !self.took_its_number || self.ended.is_some() {
            return;
        }
        a_pass_of_the_run_ended(
            self.run,
            PassThatEnded {
                pass: self.pass,
                bytes_read: self.bytes_read,
                num_bytes: self.num_bytes,
            },
        );
    }
}

/// What a read that the function of the page threw in fails with.
///
/// The kind is the one of `std::io::Error::other`: `ErrorKind::Interrupted`
/// is read again by four loops of the core, `read_line_of` of `io::vcf`,
/// `take_from` and `take_from_into` of `io::bgzf` and the `read_exact` of
/// `bytes_at` of `io::vars`, so a stop written with it would never end the
/// pass, and `ErrorKind::UnexpectedEof` is what `bytes_at` turns into the
/// error of a vars file that was cut short, so a stop written with it would
/// reach the user as a damaged file.
///
/// No user of the package reads this message: the consumer of the run throws
/// the value the function threw in place of whatever error the core made of
/// this one. It is what a run that could not be counted, the one with no
/// entry of [`RUNS`], fails with, and that run tells the page nothing, so no
/// function of an application throws in it either.
fn the_pass_was_stopped() -> std::io::Error {
    std::io::Error::other(
        "the function that is told how far a pass has got threw, and the pass ended \
         where it was reading",
    )
}

/// What a read whose bytes are more than a `u64` counts fails with, and
/// every read of that pass after it.
///
/// No target popnei builds for reaches it: a `usize` is 32 bits in wasm and
/// 64 natively, and both fit in a `u64`. What a count that silently became
/// `u64::MAX` would give a page is a bar that is full at its first read.
fn the_count_did_not_fit() -> std::io::Error {
    std::io::Error::other(
        "a read gave more bytes than the count of a pass holds, which is a defect \
         of popnei; please report it",
    )
}

impl Read for PassOverTheBytes {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        self.before_a_read()?;
        let num_read = match &mut self.bytes {
            TheBytes::InMemory(cursor) => cursor.read(buf)?,
            TheBytes::OfAFile(file) => file.read(buf)?,
        };
        self.has_read(num_read)?;
        Ok(num_read)
    }
}

impl BufRead for PassOverTheBytes {
    /// The bytes the source has ready, which for bytes in the memory of wasm
    /// are every one of them that the pass has not read, and for a file of
    /// the page what is left of the range it holds, a range being read when
    /// the pass holds none of one.
    ///
    /// What the pass has read is counted in [`BufRead::consume`], which is
    /// what the reader of a VCF takes its lines with: looking at the bytes
    /// takes none of them.
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.before_a_read()?;
        match &mut self.bytes {
            TheBytes::InMemory(cursor) => cursor.fill_buf(),
            TheBytes::OfAFile(file) => file.fill_buf(),
        }
    }

    /// `consume` gives nothing back, so the bytes it could not count end the
    /// pass at its next read instead: `has_read` keeps that in the pass.
    fn consume(&mut self, num_bytes: usize) {
        match &mut self.bytes {
            TheBytes::InMemory(cursor) => cursor.consume(num_bytes),
            TheBytes::OfAFile(file) => file.consume(num_bytes),
        }
        drop(self.has_read(num_bytes));
    }
}

impl Seek for PassOverTheBytes {
    /// Where the pass reads from next, which the reader of a vars file moves
    /// to the footer and to each batch.
    ///
    /// A seek reads no byte, so it counts nothing and tells nobody.
    fn seek(&mut self, to: SeekFrom) -> std::io::Result<u64> {
        match &mut self.bytes {
            TheBytes::InMemory(cursor) => cursor.seek(to),
            TheBytes::OfAFile(file) => file.seek(to),
        }
    }
}

/// What a source keeps in JavaScript, which nothing of Rust may hold and stay
/// `Send`.
///
/// One entry per source, which `VcfSource` and `VarsSource` hold the number
/// of, and which `free()` of the source takes out once no run over it is
/// open.
struct InJavaScript {
    /// The file the ranges of a pass are read from and the reader of them,
    /// and nothing for a source whose bytes are already in the memory of
    /// wasm.
    ///
    /// What the page holds for a source over a file is this pair and not the
    /// bytes: a `Blob` is a handle, with the name and the size of a file
    /// behind it.
    file: Option<(Blob, FileReaderSync)>,
    /// What the application is told the progress with, the function of
    /// `Variants.onProgress`, and nothing until it sets one.
    told: Option<Function>,
    /// Whether the source was freed while a run over it was open, which is
    /// when its entry goes: a pass that is reading tells the page through
    /// this entry until it is done.
    freed: bool,
}

/// One run of one consumer: which source it reads, how many passes it makes,
/// how many have begun, and what the function of the page threw.
///
/// It is taken out when the consumer returns, and for `iterBlocks` when the
/// iteration ends or the pass is freed.
struct Run {
    source: u32,
    num_passes: u32,
    passes_begun: u32,
    /// The passes of this run that are over, which the end of the run tells
    /// the page of, one call for each of them in the order of their numbers.
    ///
    /// They are put here as the reader of each pass is dropped, which is in
    /// no order of their numbers: the principal components of the variants
    /// hold both of their readers to the end and drop them together.
    passes_that_ended: Vec<PassThatEnded>,
    /// What the function that is told the progress threw, which ended a pass
    /// of this run and is what the consumer throws in place of the error the
    /// core gave.
    ///
    /// It is nothing when the run starts, so no run throws what another one
    /// was stopped with: a run that ends gives its entry back and the next
    /// one made in it is a new [`Run`].
    stopped_with: Option<JsValue>,
}

/// One pass of a run that is over: which pass of the run it was, how many
/// bytes it read and how many the file holds.
///
/// It is what the call at the end of the run carries, and it is kept in the
/// run because the pass itself is gone by then: the reader that held it was
/// dropped, which is what says that a pass is over.
#[derive(Clone, Copy)]
struct PassThatEnded {
    pass: u32,
    bytes_read: u64,
    num_bytes: u64,
}

thread_local! {
    /// What each open source keeps in JavaScript. A tab has one thread of
    /// wasm, and a table that no thread leaves is what keeps the readers of
    /// the core `Send` with a handle of JavaScript behind them.
    static IN_JAVASCRIPT: RefCell<Vec<Option<InJavaScript>>> = const { RefCell::new(Vec::new()) };
    /// The runs that are open, one per call of a consumer that has not
    /// returned.
    static RUNS: RefCell<Vec<Option<Run>>> = const { RefCell::new(Vec::new()) };
}

/// The number of the entry `what` was put in, the first free one of `table`
/// or a new one at its end, and nothing when the table holds as many entries
/// as a `u32` counts.
fn put_in<T>(table: &mut Vec<Option<T>>, what: T) -> Option<u32> {
    let at = match table.iter().position(Option::is_none) {
        Some(free) => free,
        None => {
            table.push(None);
            table.len().checked_sub(1)?
        }
    };
    // The number is taken before the entry is filled, so that a table of more
    // entries than a `u32` counts leaves the one it just made free instead of
    // holding an entry that nobody can name.
    let number = u32::try_from(at).ok()?;
    *table.get_mut(at)? = Some(what);
    Some(number)
}

/// The entry of `table` numbered `at`, and nothing when there is none.
fn entry_of<T>(table: &[Option<T>], at: u32) -> Option<&T> {
    table.get(usize::try_from(at).ok()?)?.as_ref()
}

/// The entry of `table` numbered `at`, to be changed, and nothing when there
/// is none.
fn entry_to_change<T>(table: &mut [Option<T>], at: u32) -> Option<&mut T> {
    table.get_mut(usize::try_from(at).ok()?)?.as_mut()
}

/// Where the file of a source that was just opened over the bytes of a
/// `Uint8Array` is, and the number of the entry it keeps in JavaScript, which
/// holds no file and no function to tell the progress to yet.
///
/// # Errors
///
/// When the table holds as many sources as a `u32` counts, which no tab
/// reaches: a source holds the bytes of its file in the memory of wasm.
pub(crate) fn the_bytes_of_a_new_source(
    bytes: Vec<u8>,
) -> Result<(TheFileOfASource, u32), JsPopneiError> {
    // The `Vec` wasm-bindgen filled with the bytes of the `Uint8Array` is the
    // one every pass reads: an `Arc<[u8]>` here would allocate the whole file
    // again and copy it into the new buffer, and the memory of wasm never
    // gives that back.
    let file = TheFileOfASource::InMemory(Arc::new(bytes));
    let entry = the_entry_of_a_new_source(None)?;
    Ok((file, entry))
}

/// Where the file of a source that was just opened over `file`, a file the
/// user picked in the page, is, and the number of the entry that holds it
/// with the reader of its ranges.
///
/// # Errors
///
/// When the browser has no `FileReaderSync`, which is every call outside a
/// web worker; when `Blob.size` is not a whole number of bytes popnei reads a
/// file by; and when the table holds as many sources as a `u32` counts.
pub(crate) fn the_file_of_a_new_source(
    file: Blob,
) -> Result<(TheFileOfASource, u32), JsPopneiError> {
    let reader = FileReaderSync::new().map_err(|_| {
        JsPopneiError::Refused(
            "popnei reads a `File` or a `Blob` through `FileReaderSync`, which a \
             browser gives only inside a web worker, and this call was made where \
             there is none: the main thread of a page, or node. Open the file inside \
             a web worker, or give its bytes as a `Uint8Array`."
                .to_owned(),
        )
    })?;
    let num_bytes = the_size_of_a_file(file.size())?;
    let entry = the_entry_of_a_new_source(Some((file, reader)))?;
    Ok((TheFileOfASource::OfThePage { num_bytes }, entry))
}

/// How many bytes the file of a source holds, from the `size` of its `Blob`.
///
/// `Blob.size` is a float64, so it is looked at before it becomes a count:
/// `f64::NAN as u64` is 0, which would be a file of no bytes, and a size
/// above 2^53 is a number a float64 no longer counts one by one, which the
/// two ends of every range of that file cross as.
///
/// # Errors
///
/// When the size is not a whole number of bytes from 0 to 2^53.
fn the_size_of_a_file(size: f64) -> Result<u64, JsPopneiError> {
    #[expect(
        clippy::cast_precision_loss,
        reason = "2^53 is exact as a float64, which is what it is the largest of"
    )]
    let largest = LARGEST_POSITION as f64;
    if !size.is_finite() || size < 0.0 || size > largest || size.fract() != 0.0 {
        return Err(JsPopneiError::Refused(format!(
            "this file says that it holds {size} bytes, and popnei reads a file of a \
             whole number of bytes from 0 to {LARGEST_POSITION}: the two ends of every \
             range it asks the browser for cross as a number of JavaScript, which \
             counts one by one up to that number and no further"
        )));
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "checked above to be a whole number between 0 and 2^53"
    )]
    let num_bytes = size as u64;
    Ok(num_bytes)
}

/// The number of the entry a source that was just opened keeps in JavaScript,
/// made with `file` and with no function to tell the progress to yet.
///
/// # Errors
///
/// When the table holds as many sources as a `u32` counts.
fn the_entry_of_a_new_source(file: Option<(Blob, FileReaderSync)>) -> Result<u32, JsPopneiError> {
    IN_JAVASCRIPT
        .with_borrow_mut(|sources| {
            put_in(
                sources,
                InJavaScript {
                    file,
                    told: None,
                    freed: false,
                },
            )
        })
        .ok_or_else(|| {
            JsPopneiError::Broken(
                "the page holds as many open sources of variants as a number of 32 \
                 bits counts, and this one has no entry left to be opened in"
                    .to_owned(),
            )
        })
}

/// Sets `told` as what every pass over the source numbered `source` tells the
/// page with, and takes the one that was set off when it is nothing.
pub(crate) fn tells_the_progress(source: u32, told: Option<Function>) {
    IN_JAVASCRIPT.with_borrow_mut(|sources| {
        if let Some(entry) = entry_to_change(sources, source) {
            entry.told = told;
        }
    });
}

/// The source numbered `source` was freed: its entry goes, or it is marked as
/// freed and the last run over it takes it out.
///
/// A pass that is still reading when `free()` is called reads on to its end,
/// as `docs/specs/js_sources.md` says, and it tells the page through this
/// entry while it does.
pub(crate) fn the_source_was_freed(source: u32) {
    if a_run_reads(source) {
        IN_JAVASCRIPT.with_borrow_mut(|sources| {
            if let Some(entry) = entry_to_change(sources, source) {
                entry.freed = true;
            }
        });
        return;
    }
    the_entry_of_the_source_goes(source);
}

/// Takes the entry of the source numbered `source` out of [`IN_JAVASCRIPT`],
/// which gives back the function the page set and the handle of the file the
/// ranges were read from, so that the page can let go of the file.
fn the_entry_of_the_source_goes(source: u32) {
    IN_JAVASCRIPT.with_borrow_mut(|sources| {
        if let Some(entry) = usize::try_from(source)
            .ok()
            .and_then(|at| sources.get_mut(at))
        {
            *entry = None;
        }
    });
}

/// Whether a run over the source numbered `source` is open.
fn a_run_reads(source: u32) -> bool {
    RUNS.with_borrow(|runs| runs.iter().flatten().any(|run| run.source == source))
}

/// The run of `consumer` over the source numbered `source`, which is taken
/// out of [`RUNS`] when it is dropped.
pub(crate) fn starts_a_run_of(source: u32, consumer: &Consumer) -> RunOfAConsumer {
    let run = RUNS.with_borrow_mut(|runs| {
        put_in(
            runs,
            Run {
                source,
                num_passes: consumer.num_passes(),
                passes_begun: 0,
                passes_that_ended: Vec::new(),
                stopped_with: None,
            },
        )
    });
    // A run the table had no entry left for tells the page nothing, which is
    // the whole of what it would have done: no tab reaches that many open
    // runs, and a consumer that could not be counted still gives its result.
    RunOfAConsumer(run.unwrap_or(NO_ENTRY))
}

/// The number of the pass of the run numbered `run` that is starting, 1 for
/// the first, and 0 when there is no such run.
fn the_pass_that_starts(run: u32) -> u32 {
    RUNS.with_borrow_mut(|runs| {
        let Some(run) = entry_to_change(runs, run) else {
            return 0;
        };
        run.passes_begun = run.passes_begun.saturating_add(1);
        run.passes_begun
    })
}

/// A pass of the run numbered `run` is over, which the end of that run tells
/// the page of.
///
/// A pass whose run is not in the table of this thread is not kept: its run
/// is over already, or the pass was moved to another thread, and either way
/// nobody is told of it.
fn a_pass_of_the_run_ended(run: u32, ended: PassThatEnded) {
    RUNS.with_borrow_mut(|runs| {
        if let Some(run) = entry_to_change(runs, run) {
            run.passes_that_ended.push(ended);
        }
    });
}

/// Tells the page how far each pass of `run` got, now that the run is over:
/// one call for each pass that ended, in the order of their numbers.
///
/// It is what says that a pass is over, because no read does: a pass over a
/// vars file stops after its last batch, with up to a range of bytes read
/// since the last call, and a run that fails stops wherever it failed.
/// Without it a pass over a vars file of 12231602 bytes was last told at
/// 8476400, two thirds of the way, and the bar of a page stood there.
///
/// A value thrown in one of these calls is dropped and stops nothing: the
/// run is over and there is no read left for it to end. The run is out of
/// [`RUNS`] before the first of them, so an application that starts a
/// consumer from inside one finds the table as any other call does.
fn the_run_ended(run: &Run) {
    let told = IN_JAVASCRIPT.with_borrow(|sources| entry_of(sources, run.source)?.told.clone());
    let Some(told) = told else {
        return;
    };
    let mut passes = run.passes_that_ended.clone();
    passes.sort_by_key(|ended| ended.pass);
    for ended in passes {
        // The four numbers of the `Progress` of `docs/specs/js_sources.md`,
        // as a read that tells the page sends them, with the bytes that pass
        // read where a reading pass puts the bytes it has read so far.
        let progress = Array::of4(
            &JsValue::from_f64(ended.bytes_read as f64),
            &JsValue::from_f64(ended.num_bytes as f64),
            &JsValue::from_f64(f64::from(ended.pass)),
            &JsValue::from_f64(f64::from(run.num_passes)),
        );
        drop(told.apply(&JsValue::NULL, &progress));
    }
}

/// The function of the page threw `thrown` in a pass of the run numbered
/// `run`, which is what the consumer of that run throws in place of the error
/// its read failed with.
///
/// A run whose function threw twice keeps the first value: the pass the first
/// throw was in reads no more, and a second pass of the same run is not
/// started, because the consumer gets the error of the first.
fn the_run_was_stopped(run: u32, thrown: JsValue) {
    RUNS.with_borrow_mut(|runs| {
        if let Some(run) = entry_to_change(runs, run)
            && run.stopped_with.is_none()
        {
            run.stopped_with = Some(thrown);
        }
    });
}

/// What the function of the page threw in a pass of the run numbered `run`,
/// and nothing when no pass of it was stopped.
///
/// The value is cloned out of the table, which is a handle of JavaScript
/// copied: the run keeps what stopped it for as long as it is open, and the
/// `iterBlocks` whose iteration threw it is asked for its blocks again
/// without being told of another stop.
fn what_a_run_was_stopped_with(run: u32) -> Option<JsValue> {
    RUNS.with_borrow(|runs| entry_of(runs, run)?.stopped_with.clone())
}

/// The function the page is told the progress of the run numbered `run` with,
/// and how many passes that run makes, or nothing when the run is over or its
/// source was given no function.
///
/// The function is cloned out of the table, which is a handle of JavaScript
/// copied, so that no table is borrowed while it runs.
fn what_tells_the_page(run: u32) -> Option<(Function, u32)> {
    let (source, num_passes) = RUNS.with_borrow(|runs| {
        let run = entry_of(runs, run)?;
        Some((run.source, run.num_passes))
    })?;
    let told = IN_JAVASCRIPT.with_borrow(|sources| entry_of(sources, source)?.told.clone())?;
    Some((told, num_passes))
}

/// The run a consumer holds: the number of its entry of [`RUNS`], which the
/// readers of its passes carry, and which is taken out of that table when the
/// consumer is done with it.
pub(crate) struct RunOfAConsumer(u32);

impl RunOfAConsumer {
    /// `result` as the consumer of this run gives it: an error swapped for
    /// the value the function of the page threw, when that is what ended a
    /// pass of the run.
    ///
    /// The swap is made whatever error the core gave back, so nothing
    /// depends on which reader turned the failed read into which error: the
    /// same stop is a wrong line of a VCF to one reader, a gzip member that
    /// is not there to another and a vars file that was cut short to a
    /// third, and the application is given its own value in all three.
    ///
    /// It is called while the run is open, because the value it reads goes
    /// out of [`RUNS`] with the run.
    pub(crate) fn what_the_consumer_gives<T>(
        &self,
        result: Result<T, JsPopneiError>,
    ) -> Result<T, JsPopneiError> {
        match result {
            Ok(given) => Ok(given),
            Err(error) => Err(match what_a_run_was_stopped_with(self.0) {
                Some(thrown) => JsPopneiError::Stopped(thrown),
                None => error,
            }),
        }
    }
}

impl Drop for RunOfAConsumer {
    /// Takes the run out of [`RUNS`], tells the page how far each of its
    /// passes got, and takes out the entry of a source that was freed while
    /// this was the last run reading it.
    ///
    /// Every reader of the run is dropped before this, which is what put the
    /// passes in it: a consumer drops its readers when it has its result,
    /// and the `Blocks` of an iteration holds its reader before its run and
    /// so drops it first.
    fn drop(&mut self) {
        let ended = RUNS.with_borrow_mut(|runs| {
            let at = usize::try_from(self.0).ok()?;
            runs.get_mut(at)?.take()
        });
        let Some(ended) = ended else {
            return;
        };
        the_run_ended(&ended);
        let source = ended.source;
        if a_run_reads(source) {
            return;
        }
        let freed = IN_JAVASCRIPT
            .with_borrow(|sources| entry_of(sources, source).is_some_and(|entry| entry.freed));
        if freed {
            the_entry_of_the_source_goes(source);
        }
    }
}

/// What the consumer `consumer` of `source` gives, over a run of its own.
///
/// `reads_the_source` opens the readers of the passes of the run and makes
/// the calculation over them, and the run is open while it runs and is taken
/// out of [`RUNS`] when it is over.
///
/// Every consumer but the iteration of `iterBlocks` opens its run here, so
/// that the value an application threw to stop it is what the consumer gives
/// back, in place of the error the read failed with, and no consumer has to
/// remember to make that swap itself. The iteration is the one that does not:
/// its run lives on in the `Blocks` after the call that made it has returned,
/// and [`Blocks::next_block`] is where its swap is made.
///
/// # Errors
///
/// Those of the consumer, and the value the function that is told the
/// progress threw when it stopped a pass of the run.
pub(crate) fn the_run_of<T>(
    source: &dyn OpenSource,
    consumer: &Consumer,
    reads_the_source: impl FnOnce(&RunOfAConsumer) -> Result<T, JsPopneiError>,
) -> Result<T, JsPopneiError> {
    let run = source.starts_a_run(consumer);
    let given = reads_the_source(&run);
    run.what_the_consumer_gives(given)
}

/// That the memory of wasm takes `num_bytes` more, asked for before a
/// `Uint8Array` of that length is copied into it.
///
/// The code wasm-bindgen generates for an argument of bytes asks for the
/// whole length before any code of popnei runs, and an allocation that
/// fails in wasm aborts, which is a trap: the call ends where it is and the
/// module cannot be called again, where section 11 of
/// `docs/architecture.md` asks for an `Error`. So the package calls this
/// first, with the length of the array. The memory this grew is not given
/// back to the system, which wasm cannot do, so the copy that follows finds
/// it.
///
/// # Errors
///
/// When the memory of wasm cannot take that many bytes more: a wasm module
/// addresses 4 GB, and what is already open in the tab is in those 4 GB.
#[wasm_bindgen]
pub fn room_for_bytes(num_bytes: f64) -> Result<(), JsPopneiError> {
    let no_room = || {
        JsPopneiError::NoMemory(format!(
            "the {num_bytes} bytes of this file do not fit in the memory popnei has \
             left: a page holds at most 4 GB of them at a time, and every file that \
             is open counts. A file this large is read by a program outside the \
             browser, popnei in Python among them."
        ))
    };
    // The length of a `Uint8Array` is a whole number that is not negative,
    // and everything else, a NaN among it, is refused with the same message
    // instead of being cast.
    if !num_bytes.is_finite() || num_bytes < 0.0 || num_bytes > LARGEST_ALLOCATION {
        return Err(no_room());
    }
    #[expect(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "checked above to be a number between 0 and usize::MAX"
    )]
    let wanted = num_bytes as usize;
    let mut room: Vec<u8> = Vec::new();
    room.try_reserve_exact(wanted).map_err(|_| no_room())?;
    drop(room);
    Ok(())
}

/// The most bytes one allocation of this build can hold, `usize::MAX`,
/// which in wasm is 2^32 - 1.
#[expect(
    clippy::cast_precision_loss,
    reason = "2^32 - 1 is below 2^53 and is exact as a float64"
)]
const LARGEST_ALLOCATION: f64 = usize::MAX as f64;

/// One pass over `source`, through the steps of `steps`, whose blocks hold
/// `fields` besides the genotypes, `num_vars_per_block` variants each.
///
/// The steps are taken as they are here, when the pass starts: one added
/// while it runs holds from the next pass, as `docs/specs/filters.md` says.
/// A pass with no step reads the source as it is, and it is asked for with
/// an empty `Steps` and not by leaving the argument out: a caller that could
/// omit it would read the variants of a filtered `Variants` unfiltered and
/// say nothing.
///
/// # Errors
///
/// When a name of `fields` is not a field of a block, when
/// `num_vars_per_block` is 0, when the source cannot be read again, and when
/// the genotypes of one block are more than the memory of wasm addresses.
pub(crate) fn blocks_of(
    source: &dyn OpenSource,
    fields: Vec<String>,
    num_vars_per_block: Option<usize>,
    steps: Steps,
) -> Result<Blocks, JsPopneiError> {
    let needs = needs_of_the_fields(fields.iter().map(String::as_str))?;
    // The source is asked for the size the user wants, so a reader that can
    // give it has nothing for the `Reblock` over it to cut or to join and
    // every block goes through with no copy. That `Reblock` is there for the
    // sources that give another size, the vars file whose batches were
    // written at one size and a filter among them, and it is what
    // `docs/specs/block.md` puts at the end of every `iterBlocks`.
    // The iteration is a run of one pass, which the blocks hold: a user who
    // opens twelve iterations over one source at once has twelve runs, and
    // each of them is the pass 1 of 1 of its own.
    let run = source.starts_a_run(&Consumer::IterBlocks);
    // The run is not opened by `the_run_of`, as every other consumer's is,
    // because it goes into the `Blocks` and lives on after this call. So the
    // reader is built here with the run in hand and the swap of the error is
    // made on what building it gave: the first read of the pass is the header
    // of a VCF or the footer of a vars file, made when the reader is built,
    // and the function of the page can throw in it.
    let opens_the_pass = || -> Result<Box<dyn BlockReader>, JsPopneiError> {
        let reader = source.reader(&run, num_vars_per_block)?;
        // The fields are asked of the whole chain and not of the source
        // alone: a filter asks its source for what it was asked for and for
        // the genotypes, which it needs itself.
        let mut chain = chain_of(reader, steps.steps())?;
        chain.set_needs(needs.union(Needs::GTS));
        Ok(Box::new(Reblock::new(chain, num_vars_per_block)?))
    };
    let opened = opens_the_pass();
    Ok(Blocks {
        reader: run.what_the_consumer_gives(opened)?,
        run,
        finished: false,
        num_vars: 0,
    })
}

/// How many bytes one piece of a vars file that is being written holds,
/// 1 MiB.
///
/// It is what the memory of wasm grows by at a time while a file is written
/// and what one call of `next_piece` copies out, and nothing of the format
/// depends on it: the pieces are put together in JavaScript into the one
/// array the user gets. Writing a file of 18.3 MB from a vars source grew
/// that memory by 33.3 MB with pieces of this size, by 34.2 MB with pieces
/// of 128 KiB and by 36.8 MB with pieces of 4 MiB, so what is left above
/// the file is not the pieces but what the reader and arrow-rs hold while a
/// batch is written.
const BYTES_PER_PIECE: usize = 1024 * 1024;

/// A vars file that was written, held in pieces of [`BYTES_PER_PIECE`].
///
/// The whole file is in the memory of wasm when the write is over, and that
/// memory never shrinks, so what a `Vec` that grew by doubling cost a tab
/// was the bytes it had and the buffer twice that size it copied them into,
/// at every doubling, all of it kept for as long as the page lives.
/// Measured on a vars file of 18.3 MB written from a vars source, that one
/// grew the memory of wasm by 62.4 MB and these pieces grow it by 33.3 MB.
///
/// Every piece leaves the memory of wasm as it is read, which is where the
/// copy that crosses is made, and the package puts them together into the
/// `Uint8Array` the user gets, in the heap of JavaScript.
#[wasm_bindgen]
pub struct VarsFile {
    /// The pieces in the order they were written, each of them empty once
    /// it has been given to JavaScript.
    pieces: Vec<Vec<u8>>,
    num_bytes: usize,
    /// Which piece is the next one to give.
    next: usize,
    /// The counts of the pass that wrote the file, which the package gives
    /// its user beside the bytes.
    counts: PassCounts,
}

#[wasm_bindgen]
impl VarsFile {
    /// How many bytes the whole file holds.
    #[must_use]
    pub fn num_bytes(&self) -> usize {
        self.num_bytes
    }

    /// The next piece of the file, or `undefined` when it has all been
    /// given. Each of them leaves the memory of wasm as it is read.
    pub fn next_piece(&mut self) -> Option<Vec<u8>> {
        let piece = self.pieces.get_mut(self.next)?;
        self.next = self.next.saturating_add(1);
        Some(std::mem::take(piece))
    }

    /// How many variants were written, and what each filter of the pass was
    /// given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// The sink the core writes a vars file into: it takes the bytes in pieces
/// of [`BYTES_PER_PIECE`] and never copies what it has into a larger buffer.
struct PiecesOfTheFile {
    pieces: Vec<Vec<u8>>,
    num_bytes: usize,
}

impl PiecesOfTheFile {
    fn new() -> PiecesOfTheFile {
        PiecesOfTheFile {
            pieces: Vec::new(),
            num_bytes: 0,
        }
    }

    /// The piece the next bytes go into, a new one when the last is full.
    ///
    /// # Errors
    ///
    /// When the memory of wasm cannot take another piece, which is the
    /// error of a vars file that could not be written: the core wraps what
    /// a sink says, and a tab that has no memory left for the file is told
    /// so instead of trapping on a failed allocation.
    fn piece_with_room(&mut self) -> std::io::Result<&mut Vec<u8>> {
        let full = match self.pieces.last() {
            Some(piece) => piece.len() >= piece.capacity(),
            None => true,
        };
        if full {
            let mut piece: Vec<u8> = Vec::new();
            piece.try_reserve_exact(BYTES_PER_PIECE).map_err(|_| {
                std::io::Error::new(
                    ErrorKind::OutOfMemory,
                    format!(
                        "the memory of this tab does not take {BYTES_PER_PIECE} bytes \
                         more of the file, which holds {num_bytes} bytes so far",
                        num_bytes = self.num_bytes
                    ),
                )
            })?;
            self.pieces.push(piece);
        }
        self.pieces
            .last_mut()
            .ok_or_else(|| std::io::Error::other("the file has no piece to be written into"))
    }
}

impl Write for PiecesOfTheFile {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.is_empty() {
            return Ok(0);
        }
        let piece = self.piece_with_room()?;
        let room = piece.capacity().saturating_sub(piece.len());
        let taken = room.min(bytes.len());
        let Some(head) = bytes.get(..taken) else {
            return Err(std::io::Error::other(
                "the bytes of the file are fewer than what is being taken from them",
            ));
        };
        piece.extend_from_slice(head);
        self.num_bytes = self.num_bytes.checked_add(taken).ok_or_else(|| {
            std::io::Error::other("the file holds more bytes than this machine counts")
        })?;
        Ok(taken)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Every variant of `source`, through the steps of `steps`, as a vars
/// file, one batch of `num_vars_per_block` variants after another, and
/// `None` for the size popnei chooses for the individuals of the source.
///
/// The file is built in the memory of wasm, in pieces that cross one by one:
/// a `Uint8Array` that were a view into that memory would stop being valid
/// the next time it grows. A tab holds the source and the file at once, so
/// the memory it needs is the two together.
///
/// # Errors
///
/// When `num_vars_per_block` is 0, when the source cannot be read, when a
/// block of it is not one a vars file holds, and when the memory of the tab
/// does not take the file.
pub(crate) fn bytes_of_a_vars_file(
    source: &dyn OpenSource,
    num_vars_per_block: Option<usize>,
    steps: Steps,
) -> Result<VarsFile, JsPopneiError> {
    // The source is asked for the size the batches will have, as a pass is,
    // so that the `reblock` the core puts over it has nothing to cut or to
    // join. What that saves is the memory of a block: a VCF read with the
    // size popnei chooses for 1000 individuals, 5000 variants, holds 10.4 MB
    // of genotypes while it is written, and 0.3 MB when the caller asked for
    // batches of 100. A source that cannot give that size, the vars file
    // whose batches were written at another one, leaves it to the `reblock`.
    the_run_of(source, &Consumer::WriteVars, |run| {
        let reader = source.reader(run, num_vars_per_block)?;
        // The chain of the pass stays here, lent to the core, so that the
        // counts of its filters can be read when the call is over: the loop
        // over the blocks is the core's, and so is the count of the variants
        // it wrote, which no loop of this crate sees.
        let mut chain = chain_of(reader, steps.steps())?;
        let (written, num_vars) =
            popnei::io::vars::write_vars(&mut chain, PiecesOfTheFile::new(), num_vars_per_block)?;
        Ok(VarsFile {
            pieces: written.pieces,
            num_bytes: written.num_bytes,
            next: 0,
            counts: PassCounts::of(num_vars, &chain.filtering_stats()),
        })
    })
}

/// The counts of one pass on their way to JavaScript: how many variants it
/// gave, and, for each filter of its chain, its kind, how many variants it
/// was given and how many it kept.
///
/// The filters come in the order of the chain, the outermost first, which is
/// the reverse of the order of the steps: the package turns them around, as
/// "How it runs" of the counts of `docs/specs/filters.md` says. Every count,
/// the one of the pass and the two of each filter, crosses as a number of
/// JavaScript, a float64, which holds every whole number up to 2^53: a pass
/// of wasm counts the variants of a file that is in the memory of the tab,
/// which addresses 2^32 bytes.
#[wasm_bindgen]
#[derive(Clone)]
pub struct PassCounts {
    num_vars: f64,
    kinds: Vec<String>,
    vars_processed: Vec<f64>,
    vars_kept: Vec<f64>,
}

#[wasm_bindgen]
impl PassCounts {
    /// How many variants the pass gave, after its filters.
    #[must_use]
    pub fn num_vars(&self) -> f64 {
        self.num_vars
    }

    /// The kind of each filter of the chain, the outermost first:
    /// `"missing_data"`, `"maf"`, `"obs_het"` or `"ld"`.
    #[must_use]
    pub fn kinds(&self) -> Vec<String> {
        self.kinds.clone()
    }

    /// How many variants each filter of `kinds` was given.
    #[must_use]
    pub fn vars_processed(&self) -> Vec<f64> {
        self.vars_processed.clone()
    }

    /// How many of them it kept.
    #[must_use]
    pub fn vars_kept(&self) -> Vec<f64> {
        self.vars_kept.clone()
    }
}

impl PassCounts {
    /// The `num_vars` variants of a pass and the `filtering` its chain gave,
    /// as the numbers of JavaScript.
    pub(crate) fn of(num_vars: u64, filtering: &[(&'static str, FilteringStats)]) -> PassCounts {
        PassCounts {
            num_vars: num_vars as f64,
            kinds: filtering
                .iter()
                .map(|(kind, _)| (*kind).to_owned())
                .collect(),
            vars_processed: filtering
                .iter()
                .map(|(_, stats)| stats.vars_processed as f64)
                .collect(),
            vars_kept: filtering
                .iter()
                .map(|(_, stats)| stats.vars_kept as f64)
                .collect(),
        }
    }
}

/// One pass over a source of variants, which gives them block by block.
#[wasm_bindgen]
pub struct Blocks {
    reader: Box<dyn BlockReader>,
    /// The run of the iteration, which is taken out of [`RUNS`] when the pass
    /// is freed: an iteration that a user abandons holds it until the
    /// `FinalizationRegistry` of the package frees the pass.
    ///
    /// It is read at every block, to give back the value an application threw
    /// to stop the iteration in place of the error the read failed with.
    run: RunOfAConsumer,
    /// Whether the pass is over: the reader has no more blocks, or a block
    /// was lost with an error. After either there is no block.
    finished: bool,
    /// The variants of the blocks the pass has given, which is the
    /// `num_vars` a user reads in its counts. A block that was lost with an
    /// error is not among them: it never reached the user.
    num_vars: u64,
}

#[wasm_bindgen]
impl Blocks {
    /// The next block of the pass, or `undefined` when the source has no
    /// more variants.
    ///
    /// The names of the chromosomes are taken after the block was given, as
    /// `docs/specs/block.md` says: the table of the reader grows while the
    /// file is read, and a block holds numbers of it.
    ///
    /// The block is checked before its genotypes cross as one `Int8Array`
    /// of variants x individuals x ploidy: a block whose arrays are not of
    /// its size would be read one genotype at the place of another, with
    /// nothing to show it.
    ///
    /// # Errors
    ///
    /// When a variant cannot be read, when the block is not of its own size,
    /// and when a position of the block is above [`LARGEST_POSITION`]. The
    /// block that was being built is lost with the error, and every call
    /// after it gives no block: the errors of this pass happen after the
    /// block was taken from the reader, which knows nothing of them, so it
    /// is this pass that keeps the promise of `docs/specs/block.md` that
    /// the variants after a wrong one are not handed out as if nothing had
    /// happened.
    pub fn next_block(&mut self) -> Result<Option<BlockColumns>, JsPopneiError> {
        if self.finished {
            return Ok(None);
        }
        let columns = self.columns_of_the_next_block();
        if !matches!(columns, Ok(Some(_))) {
            self.finished = true;
        }
        // The iteration is the one consumer whose run is not opened by
        // `the_run_of`, so this is where the value that stopped it takes the
        // place of the error the read failed with.
        self.run.what_the_consumer_gives(columns)
    }

    /// How many variants the pass has given, and what each filter of it was
    /// given and kept, the outermost filter first.
    ///
    /// It is read while the pass runs too, and it then holds what has been
    /// read up to there.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        PassCounts::of(self.num_vars, &self.reader.filtering_stats())
    }
}

impl Blocks {
    /// The columns of the next block of the reader, or `None` when it has no
    /// more variants. What ends the pass is [`Blocks::next_block`], which
    /// calls this one.
    fn columns_of_the_next_block(&mut self) -> Result<Option<BlockColumns>, JsPopneiError> {
        let Some(block) = self.reader.next_block()? else {
            return Ok(None);
        };
        block.check()?;
        let chroms = self.reader.chroms();
        let names = |numbers: Vec<u32>| {
            numbers
                .into_iter()
                .map(|number| {
                    chroms.name(number).map(str::to_owned).ok_or_else(|| {
                        JsPopneiError::Broken(format!(
                            "the chromosome number {number} of a block is not in the \
                             table of the reader that gave it"
                        ))
                    })
                })
                .collect::<Result<Vec<String>, JsPopneiError>>()
        };
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
        } = block;
        let alleles = alleles.map(|column| alleles_of(&column)).transpose()?;
        let (alleles, num_alleles_per_var) = match alleles {
            Some((texts, counts)) => (Some(texts), Some(counts)),
            None => (None, None),
        };
        let columns = BlockColumns {
            num_vars,
            num_individuals,
            ploidy,
            gts: Some(gts),
            chrom: chrom.map(names).transpose()?,
            pos: pos.as_deref().map(positions_of).transpose()?,
            id,
            alleles,
            num_alleles_per_var,
            qual,
        };
        // The variants of a block that is going out are counted here, after
        // every column of it crossed, where the block is the user's: a block
        // that was lost with an error of the pass itself, a position above
        // the largest a number of JavaScript holds, never reached them. A
        // variant is a row of a file, so a pass of the 18446744073709551615
        // variants this count holds is more rows than any file system
        // takes: the sum cannot reach its end. The conversion cannot fail
        // either, since a `usize` is 32 bits in wasm and 64 natively, and
        // one that did would make the count of the pass that number instead
        // of saying so.
        let num_vars_of_the_block = u64::try_from(num_vars).map_err(|_| {
            JsPopneiError::Broken(format!(
                "a block of {num_vars} variants holds more of them than the count \
                 of a pass does"
            ))
        })?;
        self.num_vars = self.num_vars.saturating_add(num_vars_of_the_block);
        Ok(Some(columns))
    }
}

/// The columns of one block on their way to TypeScript.
///
/// Each of them leaves the block the first time it is asked for, and the
/// call after that gives nothing: the package reads every column once, into
/// the object a user holds, and frees the block. What a column costs is then
/// one copy, the one wasm-bindgen makes when it crosses.
#[wasm_bindgen]
pub struct BlockColumns {
    num_vars: usize,
    num_individuals: usize,
    ploidy: usize,
    gts: Option<Vec<i8>>,
    chrom: Option<Vec<String>>,
    pos: Option<Vec<f64>>,
    id: Option<Vec<String>>,
    /// The alleles of every variant one variant after another, and how many
    /// of them each variant has: a JavaScript array of arrays of strings is
    /// not one of the types wasm-bindgen carries, so the column crosses as
    /// the core holds it, flat, and the package cuts it into the alleles of
    /// each variant.
    alleles: Option<Vec<String>>,
    num_alleles_per_var: Option<Vec<u32>>,
    qual: Option<Vec<f32>>,
}

#[wasm_bindgen]
impl BlockColumns {
    /// How many variants the block holds.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals the source has, the same for every variant.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.ploidy
    }

    /// The genotypes, `num_vars` x `num_individuals` x `ploidy` alleles,
    /// variant after variant and inside a variant individual after
    /// individual. Every block holds them, so `undefined` says that they
    /// were read already, as it does for every other column.
    pub fn gts(&mut self) -> Option<Vec<i8>> {
        self.gts.take()
    }

    /// The name of the chromosome of each variant, or `undefined` when the
    /// pass was not asked for it.
    pub fn chrom(&mut self) -> Option<Vec<String>> {
        self.chrom.take()
    }

    /// The position of each variant, 1 based as in a VCF.
    pub fn pos(&mut self) -> Option<Vec<f64>> {
        self.pos.take()
    }

    /// The id of each variant, an empty text for a variant that has none.
    pub fn id(&mut self) -> Option<Vec<String>> {
        self.id.take()
    }

    /// The alleles of every variant, the reference one first, one variant
    /// after another: the first `num_alleles_per_var()[0]` of them are the
    /// alleles of the first variant, and so on.
    pub fn alleles(&mut self) -> Option<Vec<String>> {
        self.alleles.take()
    }

    /// How many alleles each variant has, the reference one counted.
    pub fn num_alleles_per_var(&mut self) -> Option<Vec<u32>> {
        self.num_alleles_per_var.take()
    }

    /// The quality of each variant, phred scaled as the QUAL of a VCF, and
    /// NaN for a variant that has none.
    pub fn qual(&mut self) -> Option<Vec<f32>> {
        self.qual.take()
    }
}

/// The positions of a block as the numbers of JavaScript, which are float64
/// and hold a position of up to 2^53 exactly, as `docs/specs/block.md` says.
///
/// # Errors
///
/// When a position is above [`LARGEST_POSITION`]. A float64 would round it,
/// and the same file read from Python gives the position the source has, so
/// the two languages would disagree about where a variant is. No genome
/// comes near that number: the largest one known, over 1e11 bases in all of
/// its chromosomes together, is smaller by more than four orders of
/// magnitude.
pub(crate) fn positions_of(positions: &[u64]) -> Result<Vec<f64>, JsPopneiError> {
    positions
        .iter()
        .copied()
        .map(|pos| {
            if pos > LARGEST_POSITION {
                return Err(JsPopneiError::NotInJavaScript(format!(
                    "the position {pos} of a variant is above {LARGEST_POSITION}, the \
                     last whole number a number of JavaScript holds: the one after it \
                     would arrive as {LARGEST_POSITION} itself"
                )));
            }
            Ok(pos as f64)
        })
        .collect()
}

/// The alleles of every variant of a block, the reference one first, one
/// variant after another, and how many of them each variant has.
///
/// # Errors
///
/// When a variant has more alleles than a `u32` counts, which no source
/// gives: a variant of a VCF has as many alleles as its ALT column names.
fn alleles_of(column: &AllelesColumn) -> Result<(Vec<String>, Vec<u32>), JsPopneiError> {
    let mut texts = Vec::new();
    let mut counts = Vec::with_capacity(column.num_vars());
    for var in 0..column.num_vars() {
        let num_alleles = column.num_alleles(var);
        texts.extend((0..num_alleles).map(|allele| column.allele(var, allele).to_owned()));
        counts.push(u32::try_from(num_alleles).map_err(|_| {
            JsPopneiError::Broken(format!(
                "the variant {var} of a block has {num_alleles} alleles, more than a \
                 JavaScript array of counts holds"
            ))
        })?);
    }
    Ok((texts, counts))
}
