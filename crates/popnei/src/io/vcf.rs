//! The VCF reader: it gives the variants of a VCF, plain or gzipped, in
//! blocks.
//!
//! [`VcfReader::new`] takes any source of bytes and reads the header of the
//! VCF, so the individuals are known before a variant is, and
//! [`VcfReader::from_path`] does the same for a caller that has a path. The
//! reader is a [`BlockReader`]: each call gives the next run of variants as
//! the arrays of a [`Block`].
//!
//! A line is parsed into its own row of the block, so the lines of a batch
//! are parsed side by side on the threads of rayon, and in wasm, which has
//! no threads, one after another. The columns of the individuals are read
//! as bytes and never as text.
//!
//! The source is gzipped when its first two bytes are those of gzip,
//! whatever the name of the file. A VCF written by bgzip, which is what
//! nearly every gzipped VCF is, is many gzip members one after another, and
//! a decoder that stopped at the first one would give the header and no
//! variant and say nothing went wrong, so the reader decompresses with
//! flate2's `MultiGzDecoder`.
//!
//! `docs/specs/io_vcf.md` has the rules and where each one comes from.

use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

use crate::block::{
    AllelesColumn, Block, BlockReader, BlockSize, check_the_size_of_a_block,
    default_num_vars_per_block, size_of_the_blocks,
};
use crate::error::{Error, Result};
use crate::variant::{ChromTable, MAX_ALLELE, MISSING_ALLELE, Needs};

/// The ploidy a VCF is read with when the caller asks for no other, the
/// ploidy of a diploid organism. pyNei has no such argument and reports 2
/// for every file it reads.
pub const DEFAULT_PLOIDY: usize = 2;

/// Whether a VCF is read with the variants that failed a filter left out,
/// when the caller asks for no other. The owner decided in September 2026
/// that the FILTER column is honoured and that this is the default; pyNei
/// ignores that column and gives every variant.
pub const DEFAULT_ONLY_PASSED: bool = true;

/// The largest ploidy a reader takes. Nothing that has been sequenced comes
/// near it: the ploidies of the crops with the most copies of their genome
/// are 6 and 8. It is here so that a ploidy that came from a user cannot
/// ask the reader for a genotype of more alleles than a machine can hold: a
/// ploidy of `usize::MAX` asked for a vector of that many alleles at the
/// first genotype written as a dot.
pub const MAX_PLOIDY: usize = 255;

/// The two bytes every gzipped file starts with.
const GZIP_BYTES: [u8; 2] = [0x1f, 0x8b];

/// The empty gzip member of 28 bytes that bgzip writes at the end of a
/// file, which says that the file is whole: a source that bgzip wrote and
/// that does not end with it was cut short.
///
/// It is a member with the extra field `BC` of bgzip, no text and a CRC of
/// 0, and htslib writes these same bytes in every file it closes.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Where the flags of a gzip header are, and the flag that says that the
/// header carries an extra field, which is where bgzip writes its `BC`.
const GZIP_FLAGS: usize = 3;
const GZIP_HAS_AN_EXTRA_FIELD: u8 = 0x04;

/// Where the extra field of a gzip header starts when there is one: after
/// the ten bytes of the header and the two of its length.
const GZIP_EXTRA_FIELD: usize = 12;

/// The two bytes that name the extra field bgzip writes in every block of a
/// file, `BC`, and which say that the file was written by bgzip.
const BGZIP_EXTRA_FIELD: [u8; 2] = *b"BC";

/// How many bytes are read from the source before it is handed on: the two
/// of gzip, and enough of what comes after them for the message that says
/// that the source is not a VCF. They are given back in front of the
/// source, so reading them costs no byte of it.
const BYTES_LOOKED_AT: usize = 16;

/// What a column of a VCF holds when it has no value: the id of a variant
/// with no id, an ALT with no alternative allele, an allele that was not
/// called.
const MISSING_VALUE: &str = ".";

/// How many lines the reader takes from the source before it parses them,
/// on the threads of rayon.
///
/// The lines of one batch are what the threads share, and a batch ends
/// where its block does, so the last batch of a block holds what is left of
/// it and is the one the threads share worst: with 1024 lines the block of
/// 2500 variants that this reader gives for a file of 1000 diploid
/// individuals, [`crate::block::GENOTYPES_PER_BLOCK`] over its genotypes per
/// variant, is read in batches of 1024, 1024 and 452 lines, and 452 lines
/// are 25 to a thread on 18 cores. A batch of 4096 lines covers such a
/// block whole when [`BYTES_PER_BATCH`] lets it.
///
/// Measured on 21 September 2026 on the owner's Apple M5 Pro, 18 cores,
/// release, the file in the page cache, the genotypes alone and the default
/// options, the median of 5 runs of `benches/read_vcf.rs` on the VCF of
/// "Speed" of `docs/specs/io_vcf.md`, 100000 variants of 1000 individuals
/// in 403 MB, three sets of runs for each size, interleaved: 256, 1024,
/// 2048 and 4096 lines read it in 0.139, 0.105, 0.098 and 0.098 s on 18
/// threads, and bgzipped, 38 MB, where the decompression is one thread's
/// work whatever the others do, 256, 1024 and 4096 in 0.476, 0.425 and
/// 0.411 s. On one thread the sizes are the same read, plain 0.58 to 0.62 s
/// and bgzipped 0.90 to 0.95 s, which is the spread of the machine between
/// one set of runs and the next. `docs/reports/vcf-to-blocks.md` has the
/// numbers of the reader before this one, which filled one variant at a
/// time.
///
/// 4096 and 2048 are the same read because the bound in bytes is what
/// decides at both: 2500 lines of that file are 10.1 MB. 4096 is the
/// number, so that the bound in bytes is the one thing that cuts a block
/// into batches, and the benchmark takes `--lines-per-batch` for whoever
/// measures again.
#[cfg(not(target_family = "wasm"))]
const LINES_PER_BATCH: usize = 4096;

/// How many lines the reader takes from the source before it parses them
/// in wasm, where there are no threads: one, which is the reader that
/// parses a line and hands its variant out before it reads the next. A
/// batch there would hold the memory of its lines and buy nothing, since
/// the same one thread parses them.
#[cfg(target_family = "wasm")]
const LINES_PER_BATCH: usize = 1;

/// How many bytes of text a batch holds at most, whatever the number of
/// lines it was allowed.
///
/// A line of a VCF carries one genotype for every individual, so a bound
/// in lines alone lets the memory of a reader grow with the panel: the
/// review of work package 5 in `docs/reports/vcf-to-blocks.md` measured, on
/// the reader before this one and with a bound of 1024 lines and none in
/// bytes, 13.3 MB of resident memory for 1000 individuals, 82.7 MB for
/// 10000 and near 0.8 GB for 100000. The text of the lines is what that
/// memory is made of, the variants parsed from them and the growth of
/// the buffers by doubling, so bounding the text bounds all of it.
///
/// Where it cuts, a block is read in more than one batch, and the last
/// batch of a block is the one the threads of rayon share worst. Of the VCF
/// of "Speed" of `docs/specs/io_vcf.md`, 100000 variants of 1000
/// individuals in 403 MB and about 4 KB a line, the block of 2500 variants
/// is 10.1 MB of text: with 8 MiB a batch got about 2077 of those lines and
/// a block was read in two, and with 16 MiB it is read in one. Measured on
/// 21 September 2026 on the owner's Apple M5 Pro, 18 cores, release, the
/// file in the page cache, the genotypes alone and the default options, the
/// median of 5 runs of `benches/read_vcf.rs`, three sets of runs for each
/// bound, interleaved: on 18 threads 8 MiB reads the plain file in 0.098 s
/// and 16 MiB in 0.094, and the bgzipped one, 38 MB, in 0.411 s and
/// 0.392 s. On one thread the two are the same read, plain 0.58 to 0.62 s
/// and bgzipped 0.90 to 0.95 s, which is the spread of the machine between
/// one set of runs and the next.
///
/// What the 16 MiB costs is memory. The most bytes a reader held at once on
/// 18 threads, counted by an allocator that adds every allocation and
/// subtracts every free, were 23.3 MB with 8 MiB and 35.7 MB with 16 MiB on
/// that file, and 22.5 MB and 34.6 MB on a VCF of 10000 individuals, 3000
/// variants in 120 MB; the maximum resident set size of the benchmark that
/// reads it, which holds the binary and the pages the allocator keeps too,
/// went from 20.0 MB to 33.5 MB and from 24.4 MB to 32.8 MB. What is alive
/// at once is the genotypes of the block, 5.0 MB at 2500 variants of 1000
/// diploid individuals, the text of a batch, and the buffer of the file;
/// each of those grows by doubling, and a growth holds the old buffer and
/// the new one at once, so a reader of that file with one line in a batch,
/// whose text is 4 KB, already peaks at 10.2 MB.
///
/// A batch holds one line whatever its bytes are.
const BYTES_PER_BATCH: usize = 16 * 1024 * 1024;

/// How many bytes of the file [`VcfReader::from_path`] holds between two
/// calls to the file system, for the callers that give a path and not a
/// source of their own.
///
/// Every line the reader takes comes out of this buffer, and a buffer that
/// runs out is a call to the file system that the thread reading the lines
/// waits for. On 18 threads that wait is a large part of the time, since
/// the lines are read serially while the rows are parsed side by side.
/// Measured on 21 September 2026 on the owner's Apple M5 Pro, 18 cores,
/// release, the file in the page cache, the genotypes alone and the default
/// options, the median of 5 runs of `benches/read_vcf.rs` on the VCF of
/// "Speed" of `docs/specs/io_vcf.md`, 100000 variants of 1000 individuals
/// in 403 MB, with the batches this file has, [`LINES_PER_BATCH`] lines and
/// [`BYTES_PER_BATCH`]: 8 KiB reads it in 0.105 s on 18 threads, 64 KiB in
/// 0.095, 256 KiB in 0.093 and 1 MiB in 0.093. With the 1024 lines and the
/// 8 MiB the reader had before those two were measured, the same four read
/// it in 0.118, 0.108, 0.106 and 0.105 s. On one thread the four are the
/// same read, 0.58 to 0.62 s, which is the spread of the machine between
/// one set of runs and the next.
///
/// 256 KiB is where the gain stops: 1 MiB buys nothing more and holds four
/// times the bytes. `std::io::BufReader::new` would give 8 KiB.
const BYTES_OF_THE_FILE_BUFFER: usize = 256 * 1024;

/// The nine first columns of the `#CHROM` line of a VCF with genotypes. The
/// columns after them are the individuals.
const FIRST_COLUMNS: [&str; 9] = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
];

/// Where in a data line of a VCF something is wrong, which is what the
/// error of a data line carries beside the number of the line.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum VcfPlace {
    /// One of the nine first columns, by the name the `#CHROM` line gives
    /// it, `POS`.
    Column(&'static str),
    /// The column of one individual, by its name.
    Individual(String),
    /// The line as a whole, when no one column is at fault.
    Line,
}

impl fmt::Display for VcfPlace {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VcfPlace::Column(name) => write!(formatter, "the column {name}"),
            VcfPlace::Individual(name) => write!(formatter, "the column of {name}"),
            VcfPlace::Line => formatter.write_str("the line"),
        }
    }
}

/// What the caller says about the VCF it opens: the ploidy every genotype
/// of the file has, whether the variants that failed a filter are left out,
/// and how many variants a block holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfOptions {
    /// How many alleles every genotype has. 1 or more.
    pub ploidy: usize,
    /// Skip the variants whose FILTER is neither PASS nor a dot.
    pub only_passed: bool,
    /// How many variants a block holds, 1 or more, or `None` for
    /// [`default_num_vars_per_block`] for the individuals of the header.
    pub num_vars_per_block: Option<usize>,
}

impl Default for VcfOptions {
    /// [`DEFAULT_PLOIDY`], [`DEFAULT_ONLY_PASSED`] and the size of a block
    /// that popnei chooses: 2 alleles in a genotype, the variants that
    /// passed their filters alone, and as many variants in a block as the
    /// individuals of the file make room for.
    fn default() -> VcfOptions {
        VcfOptions {
            ploidy: DEFAULT_PLOIDY,
            only_passed: DEFAULT_ONLY_PASSED,
            num_vars_per_block: None,
        }
    }
}

/// A source with the bytes that were read from it to find the gzip put
/// back in front of it, so that nothing of it was consumed.
///
/// One look at the buffer of a source may give fewer bytes than the two of
/// gzip: a pipe, and a `BufRead` of a small buffer, give one at a time, and
/// looking again gives the same one, since looking does not consume. So the
/// two bytes are read and handed back here. After them the bytes come from
/// the buffer of the source itself and are not copied.
struct WithFirstBytes<R: BufRead> {
    /// The bytes read from the source before it was handed over.
    first: Vec<u8>,
    /// How many of them have been consumed.
    consumed: usize,
    source: R,
}

impl<R: BufRead> WithFirstBytes<R> {
    /// It reads `wanted` bytes of `source`, or every byte of it when it
    /// holds fewer, and gives them back in front of it.
    fn new(source: R, wanted: usize) -> std::io::Result<WithFirstBytes<R>> {
        let mut source = source;
        let mut first = Vec::with_capacity(wanted);
        while let Some(missing) = wanted.checked_sub(first.len()).filter(|left| *left > 0) {
            let buffer = source.fill_buf()?;
            if buffer.is_empty() {
                break;
            }
            let take = buffer.len().min(missing);
            first.extend_from_slice(buffer.get(..take).unwrap_or_default());
            source.consume(take);
        }
        Ok(WithFirstBytes {
            first,
            consumed: 0,
            source,
        })
    }
}

impl<R: BufRead> std::io::Read for WithFirstBytes<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let read = {
            let mut buffer = self.fill_buf()?;
            std::io::Read::read(&mut buffer, out)?
        };
        self.consume(read);
        Ok(read)
    }
}

impl<R: BufRead> BufRead for WithFirstBytes<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        let left = self.first.get(self.consumed..).unwrap_or_default();
        if left.is_empty() {
            return self.source.fill_buf();
        }
        Ok(left)
    }

    fn consume(&mut self, amount: usize) {
        if self.consumed < self.first.len() {
            self.consumed = self.consumed.saturating_add(amount).min(self.first.len());
        } else {
            self.source.consume(amount);
        }
    }
}

/// A source with the last [`BGZF_EOF`] bytes that went by to the decoder
/// kept, which is what says whether a file that bgzip wrote ends with the
/// mark of its end.
///
/// flate2 does not tell its caller where a gzip member ends, and the mark is
/// an empty member, so what a reader can see of it is the bytes themselves.
/// The source is read once and forward, so they are kept as they are
/// consumed: 28 bytes of memory whatever the file is.
struct WithTheLastBytes<R: BufRead> {
    source: R,
    /// The last bytes that were consumed, the oldest first, and fewer than
    /// 28 of them until 28 have gone by.
    last: Vec<u8>,
}

impl<R: BufRead> WithTheLastBytes<R> {
    fn new(source: R) -> WithTheLastBytes<R> {
        WithTheLastBytes {
            source,
            last: Vec::with_capacity(BGZF_EOF.len()),
        }
    }

    /// Whether the last bytes that went by are the empty block that bgzip
    /// writes at the end of a file.
    fn ends_with_the_bgzip_mark(&self) -> bool {
        self.last == BGZF_EOF
    }
}

/// `bytes` are about to be consumed: `last` keeps the last 28 of what has
/// gone by, counting them, and nothing is allocated after the first 28.
fn keep_the_last(last: &mut Vec<u8>, bytes: &[u8]) {
    let kept = bytes.len().min(BGZF_EOF.len());
    // What goes of what was kept before, which is what has no place left
    // once these bytes have theirs.
    let goes = last
        .len()
        .saturating_add(kept)
        .saturating_sub(BGZF_EOF.len());
    last.drain(..goes.min(last.len()));
    last.extend_from_slice(
        bytes
            .get(bytes.len().saturating_sub(kept)..)
            .unwrap_or_default(),
    );
}

impl<R: BufRead> std::io::Read for WithTheLastBytes<R> {
    fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
        let read = {
            let mut buffer = self.fill_buf()?;
            std::io::Read::read(&mut buffer, out)?
        };
        self.consume(read);
        Ok(read)
    }
}

impl<R: BufRead> BufRead for WithTheLastBytes<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.source.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        // The bytes that are about to go are the first `amount` of the
        // buffer, and looking at the buffer again does not consume it. A
        // buffer that cannot be filled again is one whose bytes are gone
        // with an error of the input, which the reader gives its caller.
        let WithTheLastBytes { source, last } = self;
        if let Ok(buffer) = source.fill_buf() {
            keep_the_last(
                last,
                buffer.get(..amount.min(buffer.len())).unwrap_or_default(),
            );
        }
        source.consume(amount);
    }
}

/// The bytes of the VCF, decompressed when the source was gzipped.
enum VcfSource<R: BufRead> {
    /// The source as it was given, with its first bytes in front of it.
    Plain(WithFirstBytes<R>),
    /// The source through the decoder of gzip, which goes on to the next
    /// member of the file when one ends, with the last bytes that went by
    /// to it kept and whether bgzip wrote it, which its first member says.
    ///
    /// The decoder with its buffers is 288 bytes over a `Cursor<Vec<u8>>`,
    /// where the plain source over the same bytes is 64, and a reader of a
    /// file that is not gzipped would carry the larger of the two, so this
    /// one is behind a pointer: one allocation when a gzipped file is
    /// opened, and the buffer of the lines is read through one indirection
    /// more.
    Gzipped {
        source: Box<BufReader<MultiGzDecoder<WithTheLastBytes<WithFirstBytes<R>>>>>,
        /// Whether bgzip wrote it, which is what makes the mark of the end
        /// of a bgzipped file something to ask for.
        written_by_bgzip: bool,
    },
}

impl<R: BufRead> VcfSource<R> {
    /// The bytes of `source` as the reader reads them: through the decoder
    /// of gzip when they are gzipped, and as they are when they are not.
    ///
    /// The first bytes of the source say which of the two it is and whether
    /// bgzip wrote it, and they are handed back in front of it, so nothing
    /// of it is consumed.
    ///
    /// # Errors
    ///
    /// When the source is not a VCF, which its first byte after the
    /// decompression says, and when its bytes cannot be read.
    fn of(source: R) -> Result<VcfSource<R>> {
        let source = WithFirstBytes::new(source, BYTES_LOOKED_AT)?;
        let gzipped = source.first.starts_with(&GZIP_BYTES);
        // Whether bgzip wrote the file is read from its first gzip member,
        // which carries the extra field `BC` that bgzip writes in every
        // block of a file.
        let written_by_bgzip = gzipped && written_by_bgzip(&source.first);
        let mut source = if gzipped {
            VcfSource::Gzipped {
                source: Box::new(BufReader::new(MultiGzDecoder::new(WithTheLastBytes::new(
                    source,
                )))),
                written_by_bgzip,
            }
        } else {
            VcfSource::Plain(source)
        };
        let first_bytes = source.first_bytes()?;
        if !first_bytes.starts_with(b"#") {
            return Err(Error::NotAVcf {
                found: as_text(first_bytes),
            });
        }
        Ok(source)
    }

    /// Whether bgzip wrote the source, which is what makes its bytes
    /// running out a file that was cut short and not an error of the input.
    fn written_by_bgzip(&self) -> bool {
        match self {
            VcfSource::Plain(_) => false,
            VcfSource::Gzipped {
                written_by_bgzip, ..
            } => *written_by_bgzip,
        }
    }

    /// Whether the source ended without the mark of the end of a file that
    /// bgzip wrote, which says that it was cut short. Nothing is missing
    /// from a source that bgzip did not write: it has no such mark.
    fn was_cut_short(&self) -> bool {
        match self {
            VcfSource::Plain(_) => false,
            VcfSource::Gzipped {
                source,
                written_by_bgzip,
            } => *written_by_bgzip && !source.get_ref().get_ref().ends_with_the_bgzip_mark(),
        }
    }

    /// The bytes that are already in the buffer, without consuming them.
    fn first_bytes(&mut self) -> std::io::Result<&[u8]> {
        match self {
            VcfSource::Plain(source) => source.fill_buf(),
            VcfSource::Gzipped { source, .. } => source.fill_buf(),
        }
    }

    /// The bytes of the next line, with its end of line, appended to
    /// `line`, and how many they were: 0 at the end of the source.
    ///
    /// The line is bytes and not text, because a VCF is read as bytes: the
    /// columns of the individuals are never turned into text, and whether
    /// the line is text at all is the parse of that line's to say, so that
    /// a line that is not gives the error of a data line with its number.
    fn read_line(&mut self, line: &mut Vec<u8>) -> std::io::Result<usize> {
        match self {
            VcfSource::Plain(source) => source.read_until(b'\n', line),
            VcfSource::Gzipped { source, .. } => source.read_until(b'\n', line),
        }
    }
}

/// One line of a batch that gets a row of the block, with the row it was
/// parsed into and what its parse gave.
///
/// Every line of a batch has its own, so that the threads that parse them
/// share nothing, and the buffers of the row are written over by the line
/// that is parsed into it in the next batch.
struct BatchRow {
    /// Where the bytes of the line are in the text of the batch, without
    /// its end of line.
    line: Range<usize>,
    /// Its number in the file, counted from 1 with the lines of the header.
    number: u64,
    /// What the line gave, but for its genotypes, which went into the row
    /// of the block.
    row: ParsedRow,
    /// The error of its parse, when it has one. The reader gives the error
    /// of the first line of the file that has one.
    error: Option<Error>,
}

impl BatchRow {
    /// A line with no text and no row, which the reader adds to the batch
    /// when it reads more lines than the batch has held so far.
    fn new() -> BatchRow {
        BatchRow {
            line: 0..0,
            number: 0,
            row: ParsedRow::default(),
            error: None,
        }
    }

    /// The line parsed into its row, with its genotypes into `gts`, the
    /// alleles of that variant in the block, and what went wrong into its
    /// own `error`.
    fn parse(&mut self, text: &[u8], gts: &mut [i8], rules: &RowRules<'_>) {
        #[cfg(test)]
        tests::panic_if_the_test_asked_for_it(self.number, rules);
        let line = text.get(self.line.clone()).unwrap_or_default();
        self.error = parse_row(line, self.number, rules, gts, &mut self.row).err();
    }
}

/// The lines of a batch, each parsed into the row of the block that was
/// kept for it.
///
/// Natively they are parsed on threads of rayon: no two lines write the
/// same row and none of them reads another's, and what depends on the order
/// of the file, the number of a chromosome and the place of a text in the
/// column it is appended to, is done afterwards, serially, so that neither
/// the blocks nor the numbers depend on how many threads there are. In wasm
/// there are none and the same lines are parsed one after another.
///
/// The threads are those of the pool the caller is running in, and rayon's
/// global pool, one thread per core, only when the caller is in none. That
/// is what lets a test and the benchmark read the same file on a pool of
/// one thread and on a pool of many, with `install`.
///
/// `gts` holds the rows of these lines and no other, `gts_per_variant`
/// alleles for each of them, and it is empty when the genotypes were not
/// asked for.
#[cfg(not(target_family = "wasm"))]
fn parse_rows(rows: &mut [BatchRow], text: &[u8], gts: &mut [i8], rules: &RowRules<'_>) {
    use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
    use rayon::slice::ParallelSliceMut;

    if gts.is_empty() {
        rows.par_iter_mut()
            .for_each(|row| row.parse(text, &mut [], rules));
        return;
    }
    // Every individual has one allele at least and there is one individual
    // at least, since the header of a VCF with no individual is refused, so
    // the chunks are of one allele at least, which `par_chunks_mut` asks
    // for.
    let gts_per_variant = rules.gts_per_variant().max(1);
    gts.par_chunks_mut(gts_per_variant)
        .zip(rows.par_iter_mut())
        .for_each(|(gts, row)| row.parse(text, gts, rules));
}

/// The lines of a batch, each parsed into its row, one after another, which
/// is what wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn parse_rows(rows: &mut [BatchRow], text: &[u8], gts: &mut [i8], rules: &RowRules<'_>) {
    parse_rows_one_by_one(rows, text, gts, rules);
}

/// The lines of a batch parsed one after another.
///
/// It is compiled for every target and not for wasm alone, so that the
/// cargo tests, which run natively, can parse the same lines with it and
/// with the threads and compare what the two give.
#[cfg_attr(
    all(not(target_family = "wasm"), not(test)),
    expect(
        dead_code,
        reason = "in wasm it is the parse of a batch, and natively it is what the test \
                  that compares the two ways of parsing calls; outside the tests and \
                  outside wasm nothing calls it"
    )
)]
fn parse_rows_one_by_one(rows: &mut [BatchRow], text: &[u8], gts: &mut [i8], rules: &RowRules<'_>) {
    if gts.is_empty() {
        for row in rows {
            row.parse(text, &mut [], rules);
        }
        return;
    }
    let gts_per_variant = rules.gts_per_variant().max(1);
    for (gts, row) in gts.chunks_mut(gts_per_variant).zip(rows) {
        row.parse(text, gts, rules);
    }
}

/// A reader over a VCF, which gives its variants in blocks.
///
/// It holds the individuals of the file and the names of the chromosomes of
/// the variants it has given, each with its number. It is a
/// [`BlockReader`], and `docs/specs/io_vcf.md` has the rules it reads a VCF
/// by.
pub struct VcfReader<R: BufRead + Send> {
    source: VcfSource<R>,
    options: VcfOptions,
    individuals: Vec<String>,
    chroms: ChromTable,
    needs: Needs,
    /// How many variants a block holds: the size the caller asked for, or
    /// [`default_num_vars_per_block`] for the individuals of the header.
    num_vars_per_block: usize,
    /// Which of the two that size is, which the error of a block the
    /// machine has no memory for names.
    size: BlockSize,
    /// `num_individuals` x `ploidy`, the alleles of one variant, which is
    /// the row of a block that one line is parsed into.
    gts_per_variant: usize,
    /// The bytes of the lines of the batch that are given a row, one after
    /// another. The lines that are left out are taken off it again.
    text: Vec<u8>,
    /// One for each line of the batch that is given a row, with the row it
    /// was parsed into. The ones after `filled` are of the batches before,
    /// kept for their buffers.
    batch: Vec<BatchRow>,
    /// How many lines of the batch were read and are to be parsed.
    filled: usize,
    /// How many lines the reader reads before it parses them,
    /// [`LINES_PER_BATCH`], which the tests lower to read a file in several
    /// batches.
    lines_per_batch: usize,
    /// How many bytes of text those lines hold at most,
    /// [`BYTES_PER_BATCH`]: the bound that a file of many individuals
    /// reaches before the lines are counted.
    bytes_per_batch: usize,
    /// How many batches were filled, which is how the tests see that a
    /// bound cut them.
    #[cfg(test)]
    batches_filled: u64,
    /// The number of the line that was read last, counted from 1 with the
    /// lines of the header.
    line_number: u64,
    /// The error of the line that could not be read, which is given in the
    /// place of the block it would have been in.
    line_error: Option<Error>,
    /// The error the reader gives where it would have said that there are
    /// no more variants: the source was cut short. Every variant that was
    /// read is given first.
    end_error: Option<Error>,
    /// Whether the source has been read to its end or could not be read.
    source_done: bool,
    /// Whether the reader gave its last block or an error. After either,
    /// every call gives no block.
    finished: bool,
    /// Whether a parse of a batch was begun and did not come back, which is
    /// what a panic inside the parse leaves behind: the lines of that batch
    /// were never parsed, and a reader that went on would drop them without
    /// a word.
    parsing: bool,
    /// The line whose parse panics, which the test of what a reader does
    /// after a panic in its parse sets and nothing else can.
    #[cfg(test)]
    panic_at_line: Option<u64>,
}

impl<R: BufRead + Send> VcfReader<R> {
    /// The reader over `source`, the VCF, gzipped or not, whose header it
    /// reads: the individuals are known when it returns.
    ///
    /// # Errors
    ///
    /// When the ploidy of the options is 0 or above [`MAX_PLOIDY`]; when
    /// the size of a block the caller asked for is 0 or holds more
    /// genotypes than a `usize` counts, which the size popnei chooses
    /// itself is checked for when the first block is built instead; when
    /// the source is not a VCF; when its header has not the nine first
    /// columns of a VCF with genotypes or no individual after them; when
    /// two individuals have the same name; and when the source cannot be
    /// read.
    pub fn new(source: R, options: VcfOptions) -> Result<VcfReader<R>> {
        let ploidy = options.ploidy;
        if ploidy == 0 || ploidy > MAX_PLOIDY {
            return Err(Error::VcfPloidyOutOfRange {
                ploidy,
                largest: MAX_PLOIDY,
            });
        }

        let source = VcfSource::of(source)?;
        let mut reader = VcfReader {
            source,
            options,
            individuals: Vec::new(),
            chroms: ChromTable::new(),
            needs: Needs::ALL,
            num_vars_per_block: 0,
            size: BlockSize::ChosenByPopnei,
            gts_per_variant: 0,
            text: Vec::new(),
            batch: Vec::new(),
            filled: 0,
            lines_per_batch: LINES_PER_BATCH,
            bytes_per_batch: BYTES_PER_BATCH,
            #[cfg(test)]
            batches_filled: 0,
            line_number: 0,
            line_error: None,
            end_error: None,
            source_done: false,
            finished: false,
            parsing: false,
            #[cfg(test)]
            panic_at_line: None,
        };
        reader.read_header()?;
        // The individuals are known now, so the alleles of one variant and
        // the size of a block are too. A size the caller wrote is checked
        // here, before a line of the file is read; the one popnei chooses
        // is checked when the first block is built instead, so that a file
        // opened for its individuals alone is never refused for a size that
        // nobody asked for. `docs/specs/io_vcf.md` has the case.
        let num_individuals = reader.individuals.len();
        (reader.num_vars_per_block, reader.size) = match options.num_vars_per_block {
            Some(_) => size_of_the_blocks(options.num_vars_per_block, num_individuals, ploidy)?,
            None => (
                default_num_vars_per_block(num_individuals),
                BlockSize::ChosenByPopnei,
            ),
        };
        reader.gts_per_variant = num_individuals
            .checked_mul(ploidy)
            .ok_or_else(|| reader.block_too_large(reader.num_vars_per_block))?;
        Ok(reader)
    }

    /// It skips the `##` lines and takes the individuals from the `#CHROM`
    /// line, which it leaves consumed, so that the next line read is the
    /// first data line.
    ///
    /// The line it reads into is its own: the header is read once, when the
    /// reader is built, and the text of the batch is not there yet.
    fn read_header(&mut self) -> Result<()> {
        let mut line = Vec::new();
        loop {
            line.clear();
            let number = next_line_number(self.line_number);
            let read = self.source.read_line(&mut line)?;
            if read == 0 {
                return Err(Error::VcfHeader {
                    problem: "it has no #CHROM line".to_string(),
                });
            }
            self.line_number = number;
            let Ok(text) = std::str::from_utf8(without_the_bytes_of_the_line_end(&line)) else {
                return Err(Error::VcfHeader {
                    problem: format!(
                        "the bytes of its line {number} are not valid UTF-8, and a VCF is text"
                    ),
                });
            };
            if text.starts_with("##") {
                continue;
            }
            self.individuals = individuals_of(text, self.line_number)?;
            return Ok(());
        }
    }

    /// The error of a block of `num_vars_per_block` variants of this file
    /// that the machine does not give the memory for.
    fn block_too_large(&self, num_vars_per_block: usize) -> Error {
        Error::BlockTooLarge {
            num_vars_per_block,
            num_individuals: self.individuals.len(),
            ploidy: self.options.ploidy,
            size: self.size,
        }
    }
}

impl<R: BufRead + Send> VcfReader<R> {
    /// How many lines the reader takes from the source before it parses
    /// them, which is [`LINES_PER_BATCH`] until this is called. A batch
    /// holds one line at least, whatever this says, and no more lines with
    /// a row than the block it is filling has room for.
    ///
    /// It is hidden from the documentation and it is not part of what
    /// popnei promises: it is for the benchmark `benches/read_vcf.rs`,
    /// which times a file with one batch after another, and for the tests,
    /// which read a file of a few hundred lines in several batches and one
    /// line at a time, the batch of wasm. What a read gives does not depend
    /// on it.
    #[doc(hidden)]
    pub fn set_lines_per_batch(&mut self, lines: usize) {
        self.lines_per_batch = lines.max(1);
    }

    /// How many bytes of text the reader takes from the source before it
    /// parses what it read, which is [`BYTES_PER_BATCH`] until this is
    /// called. A batch holds one line at least, whatever this says.
    ///
    /// Hidden and outside what popnei promises, like
    /// [`VcfReader::set_lines_per_batch`], and for the same two callers.
    #[doc(hidden)]
    pub fn set_bytes_per_batch(&mut self, bytes: usize) {
        self.bytes_per_batch = bytes.max(1);
    }

    /// How many batches the reader has filled, which is what says whether a
    /// bound cut them.
    #[cfg(test)]
    fn batches_filled(&self) -> u64 {
        self.batches_filled
    }

    /// The line whose parse panics, for the test of what a reader does
    /// after a panic inside its parse. No VCF panics the parse.
    #[cfg(test)]
    fn panic_at_line(&mut self, line: u64) {
        self.panic_at_line = Some(line);
    }

    /// An empty column reserved for `num_items`, or `None` when `field` is
    /// not among what the blocks of this reader hold.
    ///
    /// The memory is asked for with `try_reserve_exact`, which gives it
    /// back as an error: `Vec::with_capacity` ends the process when the
    /// machine has not the memory, and panics above what a `Vec` holds, and
    /// a size that a caller of popnei wrote reaches both.
    fn reserved_column<T>(&self, field: Needs, num_items: usize) -> Result<Option<Vec<T>>> {
        if !self.needs.contains(field) {
            return Ok(None);
        }
        let mut column = Vec::new();
        column
            .try_reserve_exact(num_items)
            .map_err(|_| self.block_too_large(self.num_vars_per_block))?;
        Ok(Some(column))
    }

    /// An empty block with every column the reader was asked for, each
    /// reserved for a full block.
    ///
    /// # Errors
    ///
    /// When the genotypes of a full block are more than a `usize` counts,
    /// which is where the size that popnei chose itself is checked, and
    /// when the machine does not give the memory of one of the columns.
    fn start_block(&self) -> Result<Block> {
        let num_vars = self.num_vars_per_block;
        let too_large = || self.block_too_large(num_vars);
        // The genotypes of a full block, which is where the size that
        // popnei chose is checked: the caller asked for none, so nothing
        // was refused when the reader was built.
        let gts_per_block = check_the_size_of_a_block(
            num_vars,
            self.individuals.len(),
            self.options.ploidy,
            self.size,
        )?;
        let mut gts = Vec::new();
        if self.needs.contains(Needs::GTS) {
            gts.try_reserve_exact(gts_per_block)
                .map_err(|_| too_large())?;
        }
        let alleles = match self.needs.contains(Needs::ALLELES) {
            true => Some(AllelesColumn::with_num_vars(num_vars).map_err(|_| too_large())?),
            false => None,
        };
        Ok(Block {
            num_vars: 0,
            num_individuals: self.individuals.len(),
            ploidy: self.options.ploidy,
            gts,
            chrom: self.reserved_column(Needs::CHROM_POS, num_vars)?,
            pos: self.reserved_column(Needs::CHROM_POS, num_vars)?,
            id: self.reserved_column(Needs::ID, num_vars)?,
            alleles,
            qual: self.reserved_column(Needs::QUAL, num_vars)?,
        })
    }

    /// It reads the next lines of the source into the batch, over the text
    /// and the rows of the batch before, and keeps the ones that are given
    /// a row of the block: an empty line has none, and neither has a line
    /// whose FILTER failed when the options leave those out.
    ///
    /// It keeps `room` lines at most, so that a batch ends where its block
    /// does, and stops at the two bounds of a batch or at the end of the
    /// source. A line that cannot be read ends the batch and its error is
    /// kept apart, to be given in the place of the block it would have been
    /// in.
    fn fill_batch(&mut self, room: usize) {
        let VcfReader {
            source,
            options,
            text,
            batch,
            filled,
            lines_per_batch,
            bytes_per_batch,
            line_number,
            line_error,
            end_error,
            source_done,
            #[cfg(test)]
            batches_filled,
            ..
        } = self;
        *filled = 0;
        text.clear();
        #[cfg(test)]
        {
            *batches_filled = batches_filled.saturating_add(1);
        }
        // The text of the lines that were read, which bounds the batch
        // beside their number: one line of a file of many individuals is
        // where the memory of a reader would otherwise grow without a
        // bound.
        let mut bytes: usize = 0;
        let mut lines_read: usize = 0;
        while lines_read < *lines_per_batch && bytes < *bytes_per_batch && *filled < room {
            let number = next_line_number(*line_number);
            let start = text.len();
            match source.read_line(text) {
                Ok(0) => {
                    *source_done = true;
                    // The source is at its end, which is where the mark of
                    // the end of a file that bgzip wrote has to be: the
                    // error waits for the variants that were read to be
                    // given.
                    if source.was_cut_short() {
                        *end_error = Some(Error::VcfBgzipEndMissing);
                    }
                    break;
                }
                Ok(read) => {
                    *line_number = number;
                    // A batch holds `lines_per_batch` lines at most, so the
                    // count does not reach the largest `usize`, and the
                    // bytes of a batch stop growing at the line that
                    // reaches the bound.
                    lines_read = lines_read.saturating_add(1);
                    bytes = bytes.saturating_add(read);
                }
                Err(error) => {
                    *line_number = number;
                    *source_done = true;
                    // A source that bgzip wrote and whose bytes ran out was
                    // cut short, wherever the cut falls: a cut inside a
                    // gzip member is this error of the decoder, and one
                    // where a member ends is the mark that is not at the
                    // end of the file. A user whose download stopped is
                    // told the same thing by both. An error of the file
                    // system, a disc that fails, has another kind and stays
                    // the error of the input it is.
                    if source.written_by_bgzip() && the_bytes_ran_out(&error) {
                        *end_error = Some(Error::VcfBgzipEndMissing);
                    } else {
                        *line_error = Some(Error::Io(error));
                    }
                    break;
                }
            }
            // The serial pass that says which lines are given a row. It
            // gives no error: a line with fewer than seven columns has no
            // FILTER to find and is given one, whose parse gives the error
            // of a line with too few columns.
            let read_bytes = text.get(start..).unwrap_or_default();
            let end = start.saturating_add(without_the_bytes_of_the_line_end(read_bytes).len());
            let line = text.get(start..end).unwrap_or_default();
            if line.is_empty() || (options.only_passed && !filter_passed(line)) {
                text.truncate(start);
                continue;
            }
            if batch.len() <= *filled {
                batch.push(BatchRow::new());
            }
            let Some(row) = batch.get_mut(*filled) else {
                break;
            };
            row.line = start..end;
            row.number = number;
            row.error = None;
            *filled = filled.saturating_add(1);
        }
    }

    /// The rows of the lines of the batch appended to the block, after the
    /// variants it holds already, in the order of the file, and how many
    /// they were.
    ///
    /// The count of the variants of the block is not written here:
    /// `build_block` writes it on the one path where the block lives, and
    /// on every other the block is lost with an error.
    ///
    /// The genotypes are in the block already: the parse wrote them
    /// straight into their rows. What is appended here is the columns,
    /// which are one buffer each, and the numbers of the chromosomes, which
    /// follow the order of the variants that are given and not the order in
    /// which the lines were parsed.
    ///
    /// # Errors
    ///
    /// The error of the first line of the batch whose parse failed, which
    /// is the first one of the file, since the batches are parsed one after
    /// another. The block is lost with it.
    fn append_batch(&mut self, block: &mut Block) -> Result<usize> {
        let VcfReader {
            chroms,
            batch,
            filled,
            ..
        } = self;
        let rows = batch.get_mut(..*filled).unwrap_or_default();
        for line in rows.iter_mut() {
            if let Some(error) = line.error.take() {
                return Err(error);
            }
            let row = &line.row;
            if let Some(chrom) = block.chrom.as_mut() {
                chrom.push(chroms.intern(&row.chrom));
            }
            if let Some(pos) = block.pos.as_mut() {
                pos.push(row.pos);
            }
            if let Some(id) = block.id.as_mut() {
                id.push(row.id.clone());
            }
            if let Some(alleles) = block.alleles.as_mut() {
                alleles.push(row.alleles());
            }
            if let Some(qual) = block.qual.as_mut() {
                qual.push(row.qual);
            }
        }
        Ok(*filled)
    }

    /// The next block: batches of lines read and parsed into its rows until
    /// it holds the variants it was asked for, or until the source ends.
    ///
    /// # Errors
    ///
    /// The ones of [`VcfReader::start_block`], the error of a line that
    /// could not be read, and the error of the first wrong line of the
    /// file. The block that was being built is lost with any of them.
    fn build_block(&mut self) -> Result<Option<Block>> {
        let mut block = self.start_block()?;
        let mut num_vars: usize = 0;
        while num_vars < self.num_vars_per_block {
            if let Some(error) = self.line_error.take() {
                return Err(error);
            }
            if self.source_done {
                break;
            }
            // The room left in the block, which bounds the batch beside its
            // own two bounds: a batch ends where its block does.
            let room = self.num_vars_per_block.saturating_sub(num_vars);
            self.fill_batch(room);
            // The rows of the lines of this batch, after the ones that are
            // in the block already. They are filled with the missing allele
            // and the parse writes every one of them: a row it did not
            // write is in a block that is lost with an error.
            let start_of_the_rows = block.gts.len();
            if self.needs.contains(Needs::GTS) {
                let alleles = self
                    .filled
                    .checked_mul(self.gts_per_variant)
                    .and_then(|alleles| start_of_the_rows.checked_add(alleles))
                    .ok_or_else(|| self.block_too_large(self.num_vars_per_block))?;
                block.gts.resize(alleles, MISSING_ALLELE);
            }
            let VcfReader {
                options,
                individuals,
                needs,
                text,
                batch,
                filled,
                parsing,
                #[cfg(test)]
                panic_at_line,
                ..
            } = self;
            let rules = RowRules {
                needs: *needs,
                ploidy: options.ploidy,
                individuals,
                #[cfg(test)]
                panic_at_line: *panic_at_line,
            };
            let rows = batch.get_mut(..*filled).unwrap_or_default();
            let gts = block.gts.get_mut(start_of_the_rows..).unwrap_or_default();
            // A panic of a worker unwinds through here, and what says so
            // afterwards is this flag, which is set until the parse comes
            // back.
            *parsing = true;
            parse_rows(rows, text, gts, &rules);
            *parsing = false;
            let given = self.append_batch(&mut block)?;
            num_vars = num_vars.saturating_add(given);
        }
        if num_vars == 0 {
            // Where the reader says that there are no more variants is
            // where a source that was cut short says so, after every
            // variant it did hold was given.
            return match self.end_error.take() {
                Some(error) => Err(error),
                None => Ok(None),
            };
        }
        block.num_vars = num_vars;
        Ok(Some(block))
    }
}

/// Whether the line is given a row: its FILTER, the bytes between its sixth
/// and its seventh tab, is `PASS` or a dot, which says that no filter was
/// applied to it.
///
/// A line with fewer than seven columns has no FILTER to find and is given
/// a row, so that its parse gives the error of a line with too few columns,
/// the same one whatever the options say.
fn filter_passed(line: &[u8]) -> bool {
    let mut tabs = memchr::memchr_iter(b'\t', line);
    let Some(sixth) = tabs.nth(5) else {
        return true;
    };
    let start = sixth.saturating_add(1);
    let end = tabs.next().unwrap_or(line.len());
    let filter = line.get(start..end).unwrap_or_default();
    filter == b"PASS" || filter == MISSING_VALUE.as_bytes()
}

impl<R: BufRead + Send> fmt::Debug for VcfReader<R> {
    /// What the reader was built with and where it has got to. The source
    /// is left out, so that a reader over a source that has no `Debug` has
    /// one.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("VcfReader")
            .field("individuals", &self.individuals)
            .field("options", &self.options)
            .field("needs", &self.needs)
            .field("num_vars_per_block", &self.num_vars_per_block)
            .field("chroms", &self.chroms)
            .field("line_number", &self.line_number)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl<R: BufRead + Send> BlockReader for VcfReader<R> {
    /// The next block of the file, which holds one variant at least, and
    /// `None` when there are no more variants and at every call after that.
    ///
    /// # Errors
    ///
    /// When a line of the source cannot be read, when a data line is not
    /// one popnei reads, and when the machine does not give the memory of
    /// the block. The block that was being built is lost with the error,
    /// the blocks before it were given, and every call after it gives
    /// `None`.
    fn next_block(&mut self) -> Result<Option<Block>> {
        if self.finished {
            return Ok(None);
        }
        // A panic inside the parse of a batch unwinds through the reader
        // and leaves this: the lines of that batch were never parsed. The
        // error is given once, as the error of a reader is, and the calls
        // after it give no block.
        if self.parsing {
            self.finished = true;
            return Err(Error::VcfParseNotFinished {
                line: self.line_number,
            });
        }
        match self.build_block() {
            Ok(Some(block)) => Ok(Some(block)),
            Ok(None) => {
                self.finished = true;
                Ok(None)
            }
            Err(error) => {
                self.finished = true;
                Err(error)
            }
        }
    }

    fn individuals(&self) -> &[String] {
        &self.individuals
    }

    fn ploidy(&self) -> usize {
        self.options.ploidy
    }

    fn chroms(&self) -> &ChromTable {
        &self.chroms
    }

    /// The fields the blocks from the next one on hold. Every line of a
    /// batch is parsed into the block that asked for it, so nothing is
    /// parsed and not given when this is called.
    fn set_needs(&mut self, needs: Needs) {
        self.needs = needs;
    }
}

impl VcfReader<BufReader<File>> {
    /// The reader over the VCF at `path`, gzipped or not, for the callers
    /// that have a path and not bytes.
    ///
    /// # Errors
    ///
    /// When the file cannot be opened, with the path in the error, and
    /// everything [`VcfReader::new`] fails with.
    pub fn from_path(path: &Path, options: VcfOptions) -> Result<Self> {
        let file = File::open(path).map_err(|error| Error::FileNotOpened {
            path: path.to_path_buf(),
            source: error,
        })?;
        VcfReader::new(
            BufReader::with_capacity(BYTES_OF_THE_FILE_BUFFER, file),
            options,
        )
    }
}

/// The individuals of the `#CHROM` line, whose nine first columns have to
/// be the nine of a VCF with genotypes.
fn individuals_of(chrom_line: &str, line_number: u64) -> Result<Vec<String>> {
    let columns: Vec<&str> = chrom_line.split('\t').collect();
    for (index, expected) in FIRST_COLUMNS.iter().enumerate() {
        match columns.get(index) {
            Some(found) if found == expected => {}
            Some(found) => {
                return Err(Error::VcfHeader {
                    problem: format!(
                        "its line {line_number} has `{found}` where a VCF with genotypes \
                         has `{expected}`; the nine first columns are {first}",
                        first = FIRST_COLUMNS.join(" "),
                    ),
                });
            }
            None => {
                return Err(Error::VcfHeader {
                    problem: format!(
                        "its line {line_number} has no `{expected}` column; the nine \
                         first columns of a VCF with genotypes are {first}, and one \
                         column per individual comes after them",
                        first = FIRST_COLUMNS.join(" "),
                    ),
                });
            }
        }
    }
    let individuals = columns.get(FIRST_COLUMNS.len()..).unwrap_or_default();
    if individuals.is_empty() {
        return Err(Error::VcfHeader {
            problem: format!(
                "its line {line_number} has a FORMAT column and no individual after it"
            ),
        });
    }
    let mut seen = HashSet::with_capacity(individuals.len());
    for name in individuals {
        if name.is_empty() {
            return Err(Error::VcfHeader {
                problem: format!(
                    "an individual of its line {line_number} has no name, which is what \
                     a #CHROM line that ends in a tab has"
                ),
            });
        }
        if !seen.insert(*name) {
            return Err(Error::VcfHeader {
                problem: format!("two individuals of its line {line_number} are called `{name}`"),
            });
        }
    }
    Ok(individuals.iter().map(|name| (*name).to_string()).collect())
}

/// Whether an error of the input is the bytes of the source running out in
/// the middle of what was being read.
///
/// It is what flate2 gives for a gzip member that is cut short, "incomplete
/// deflate stream" while a member is being decompressed and "unexpected end
/// of file" in the header or the trailer of one, both of the kind
/// `UnexpectedEof`. An error of the file system has another kind, so a disc
/// that fails is not read as a file that was cut short.
fn the_bytes_ran_out(error: &std::io::Error) -> bool {
    error.kind() == std::io::ErrorKind::UnexpectedEof
}

/// Whether bgzip wrote the source, which its first gzip member says: bgzip
/// writes the extra field `BC` in the header of every block it makes, and
/// nothing else does.
///
/// `first` are the first bytes of the source, which is all that is looked
/// at: bgzip writes that field first and writes no other, so a header whose
/// extra field starts with another one is not one of bgzip's.
fn written_by_bgzip(first: &[u8]) -> bool {
    let has_an_extra_field = first
        .get(GZIP_FLAGS)
        .is_some_and(|flags| flags & GZIP_HAS_AN_EXTRA_FIELD != 0);
    let names_the_field = first
        .get(GZIP_EXTRA_FIELD..GZIP_EXTRA_FIELD.saturating_add(BGZIP_EXTRA_FIELD.len()))
        .is_some_and(|name| name == BGZIP_EXTRA_FIELD);
    has_an_extra_field && names_the_field
}

/// The number of the next line. A file of `u64::MAX` lines cannot be
/// written, so the saturation is not reached and no error is worth adding
/// for it.
fn next_line_number(line_number: u64) -> u64 {
    line_number.saturating_add(1)
}

/// Bytes of a file as a message shows them, between backticks, with the
/// ones that are not printable as their number, `\\xff`: a byte that is not
/// text shown as text is the one character that stands for everything that
/// could not be read, and a user looking for it in their file needs the
/// byte.
///
/// The first [`BYTES_LOOKED_AT`] bytes at most, which is what a column of
/// an individual, an allele number and the start of a file are worth
/// showing of.
fn shown(bytes: &[u8]) -> String {
    let mut text = String::new();
    for byte in bytes.iter().take(BYTES_LOOKED_AT) {
        if byte.is_ascii_graphic() || *byte == b' ' {
            text.push(char::from(*byte));
        } else {
            text.push_str(&format!("\\x{byte:02x}"));
        }
    }
    format!("`{text}`")
}

/// The first bytes of a source, for the message that says it is not a VCF.
fn as_text(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "nothing: it holds no byte".to_string();
    }
    shown(bytes)
}

/// The next column of a data line, or the error that says the line is not
/// one of a VCF with genotypes.
fn next_column<'a>(
    columns: &mut impl Iterator<Item = &'a str>,
    name: &str,
    line: u64,
) -> Result<&'a str> {
    columns.next().ok_or_else(|| Error::VcfDataLine {
        line,
        place: VcfPlace::Line,
        problem: format!(
            "it has no {name} column; a VCF with genotypes has the nine columns \
             {first} and one column per individual after them",
            first = FIRST_COLUMNS.join(" "),
        ),
    })
}

/// The position of the variant, 1 based as in the VCF.
fn parse_position(text: &str, line: u64) -> Result<u64> {
    text.parse().map_err(|_| Error::VcfDataLine {
        line,
        place: VcfPlace::Column("POS"),
        problem: format!("`{text}` is not a position"),
    })
}

/// The quality of the variant, `None` when the column is a dot.
///
/// A quality that is a number and not a finite one is an error of its
/// column: NaN is what a block holds for a variant with no quality, so a
/// `nan` in the file would be read as a variant that has none, and an
/// infinite quality is a probability of no variant of 0, which is not what
/// phred scaling says. `inf`, `1e400`, which is above what a float of 64
/// bits holds, and `1e39`, which is above what the 32 bits of the column of
/// a block hold, all read as an infinite one. The owner decided this on 20
/// September 2026; the option not taken was to keep what the float gave,
/// which is what pyNei does.
fn parse_quality(text: &str, line: u64) -> Result<Option<f32>> {
    if text == MISSING_VALUE {
        return Ok(None);
    }
    let wrong = |problem: String| Error::VcfDataLine {
        line,
        place: VcfPlace::Column("QUAL"),
        problem,
    };
    let Ok(quality) = text.parse::<f32>() else {
        return Err(wrong(format!("`{text}` is not a quality")));
    };
    if !quality.is_finite() {
        return Err(wrong(format!(
            "`{text}` is not a finite quality, and a quality is minus ten times the \
             base ten logarithm of the probability that there is no variant at that \
             site; a variant with no quality has a dot there"
        )));
    }
    Ok(Some(quality))
}

/// How many alleles REF and ALT declare, the reference and the alternative
/// ones.
///
/// An allele with no letter in it, which is what an empty REF, an empty ALT
/// or an ALT that ends in a comma gives, is an error of its column:
/// bcftools reads no alternative allele in `T,`, and an allele of no
/// letters would take a number that a genotype could then carry.
fn count_alleles(reference: &str, alternatives: &str, line: u64) -> Result<usize> {
    let wrong = |column: &'static str| Error::VcfDataLine {
        line,
        place: VcfPlace::Column(column),
        problem: "it has an allele with no letter in it".to_string(),
    };
    if reference.is_empty() {
        return Err(wrong("REF"));
    }
    if allele_texts(reference, alternatives).any(str::is_empty) {
        return Err(wrong("ALT"));
    }
    Ok(allele_texts(reference, alternatives).count())
}

/// The texts of the alleles of a variant, the reference first and then the
/// alternative ones, which are none when ALT is a dot.
fn allele_texts<'a>(reference: &'a str, alternatives: &'a str) -> impl Iterator<Item = &'a str> {
    let alternatives = if alternatives == MISSING_VALUE {
        None
    } else {
        Some(alternatives)
    };
    std::iter::once(reference).chain(alternatives.into_iter().flat_map(|texts| texts.split(',')))
}

/// Where `GT` is among the keys of the FORMAT column, which is where the
/// genotype of each individual is in its own column.
fn gt_index_of(format: &str, line: u64) -> Result<usize> {
    format
        .split(':')
        .position(|key| key == "GT")
        .ok_or_else(|| Error::VcfDataLine {
            line,
            place: VcfPlace::Column("FORMAT"),
            problem: format!("`{format}` has no GT key, and GT is the genotype"),
        })
}

/// What the parse of a data line into a row of a block needs to know, which
/// is the same for every line of a file.
struct RowRules<'a> {
    /// Which fields are parsed. A column that is not asked for is not read
    /// and not checked.
    needs: Needs,
    /// How many alleles every genotype of the file holds.
    ploidy: usize,
    /// The individuals of the header, in the order of their columns, by the
    /// name that an error of one of them carries.
    individuals: &'a [String],
    /// The line whose parse panics. No VCF makes the parse panic, and this
    /// is how the test of what a reader does after a panic in its parse
    /// makes one happen; nothing outside the tests can set it.
    #[cfg(test)]
    panic_at_line: Option<u64>,
}

impl RowRules<'_> {
    /// The alleles of one variant, the individuals times the ploidy, which
    /// is the row of a block that one line is parsed into.
    ///
    /// A reader refuses a file whose individuals times its ploidy are more
    /// than a `usize` counts before it reads a line, so the saturation is
    /// not reached; what does reach it is a caller of the row parser with a
    /// defect, and the parse then refuses the row it was given, whose
    /// length cannot be what it asks for.
    fn gts_per_variant(&self) -> usize {
        self.individuals.len().saturating_mul(self.ploidy)
    }
}

/// One row of a block as one data line gives it, but for the genotypes,
/// which go straight into the row of the block.
///
/// The chromosome is here as the name it has in the line: the number it
/// gets belongs to the order of the variants that are given, and the reader
/// gives it serially, after the rows of a batch were parsed side by side.
/// The id and the alleles are here and not in the block for the same
/// reason: the texts of a column are one buffer, which the rows are
/// appended to in order.
///
/// The buffers of a row are written over by the next line parsed into it,
/// so a reader that keeps one row for each line of a batch allocates
/// nothing after its first batch.
#[derive(Debug, Default)]
struct ParsedRow {
    /// The name of the chromosome, empty when the chromosome and the
    /// position were not asked for.
    chrom: String,
    /// The position, 1 based as in the VCF, and 0 when it was not asked for.
    pos: u64,
    /// The id of the variant, empty when it has none and when it was not
    /// asked for.
    id: String,
    /// The texts of the alleles, the reference first. Only the first
    /// [`ParsedRow::num_allele_texts`] of them are of this line: the strings
    /// after them are the ones of the lines parsed into this row before, and
    /// they are kept so that a line with more alleles than the line before
    /// it writes over a string instead of allocating one.
    alleles: Vec<String>,
    /// How many of the strings of `alleles` are of this line.
    num_allele_texts: usize,
    /// The quality, NaN when the variant has none and when it was not asked
    /// for.
    qual: f32,
    /// How many alleles REF and ALT declare, which is counted for every line
    /// that is parsed, whether or not the texts of the alleles are kept,
    /// because it is what says whether an allele number of a genotype is one
    /// of the alleles of the variant.
    num_alleles: usize,
}

impl ParsedRow {
    /// The texts of the alleles of the line that was parsed into this row,
    /// the reference first, and none when the alleles were not asked for.
    fn alleles(&self) -> &[String] {
        self.alleles
            .get(..self.num_allele_texts)
            .unwrap_or_default()
    }

    /// The row with nothing of the line that was parsed into it before,
    /// which the parse of a line calls before it writes anything.
    fn clear(&mut self) {
        self.chrom.clear();
        self.pos = 0;
        self.id.clear();
        self.num_allele_texts = 0;
        self.qual = f32::NAN;
        self.num_alleles = 0;
    }

    /// The texts of the alleles of REF and ALT, written over the strings the
    /// row holds.
    fn fill_alleles(&mut self, reference: &str, alternatives: &str) {
        for text in allele_texts(reference, alternatives) {
            match self.alleles.get_mut(self.num_allele_texts) {
                Some(allele) => {
                    allele.clear();
                    allele.push_str(text);
                }
                None => self.alleles.push(text.to_string()),
            }
            // The alleles of one line are at most its bytes, so the count
            // does not reach the largest `usize`.
            self.num_allele_texts = self.num_allele_texts.saturating_add(1);
        }
    }
}

/// The columns of a part of a data line, cut at the tabs, as bytes.
///
/// The tabs are searched for with `memchr`, which reads the bytes a machine
/// word at a time, and not with a loop over one byte after another, which is
/// what `[u8]::split` does. Which of the two is faster depends on how long a
/// column is, and the columns of the individuals are where nearly all the
/// bytes of a VCF with genotypes are. Both were timed on the owner's Apple
/// M5 Pro, a release build, one thread, the lines in memory and the
/// genotypes alone asked for, the median of three runs, on the two shapes a
/// column of an individual has:
///
/// | the columns of the individuals | `memchr` | a loop over the bytes |
/// |---|---|---|
/// | `0/1`, the FORMAT `GT`, 100000 lines x 1000 individuals | 0.592 s | 0.560 s |
/// | `0/1:20,30:50:99:0,120,1800`, the FORMAT `GT:AD:DP:GQ:PL`, 10000 lines x 1000 individuals | 0.069 s | 0.110 s |
///
/// The first file is `crates/popnei/benches/make_big_vcf.py`'s, the one of
/// "Speed" of `docs/specs/io_vcf.md`, and the second is its first 10000
/// lines with the four other values that a VCF of a variant caller carries
/// added to every column. `memchr` costs 6 in 100 on columns of 3 bytes and
/// gives 1.6 times on the columns of 26 that a called VCF has.
struct ByteColumns<'a> {
    /// The bytes the columns are cut from.
    bytes: &'a [u8],
    /// Where in `bytes` the column that comes next starts.
    start: usize,
    /// The tabs of `bytes` that have not been reached yet.
    tabs: memchr::Memchr<'a>,
    /// Whether the last column was given: the bytes after the last tab are
    /// one column, and after it there are none.
    done: bool,
}

impl<'a> ByteColumns<'a> {
    /// The columns of `bytes`: one when there is no tab in them, and one
    /// empty column when there are no bytes at all, since a column is what
    /// lies between two tabs and the bytes of a line begin and end one.
    fn new(bytes: &'a [u8]) -> ByteColumns<'a> {
        ByteColumns {
            bytes,
            start: 0,
            tabs: memchr::memchr_iter(b'\t', bytes),
            done: false,
        }
    }
}

impl<'a> Iterator for ByteColumns<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<&'a [u8]> {
        if self.done {
            return None;
        }
        let end = self.tabs.next().unwrap_or(self.bytes.len());
        let column = self.bytes.get(self.start..end)?;
        match end.checked_add(1) {
            Some(after) if after <= self.bytes.len() => self.start = after,
            _ => self.done = true,
        }
        Some(column)
    }
}

/// The data line `line`, the line `number` of the file, parsed into one row
/// of a block: its genotypes into `gts`, the `num_individuals` x `ploidy`
/// alleles of that row, and the rest of its fields into `row`.
///
/// Nothing is shared between two lines, so the lines of a batch are parsed
/// side by side, each into its own row. `gts` holds one allele for each
/// individual of `rules` times the ploidy, and no allele when the genotypes
/// were not asked for; the alleles of individual `i` are `gts[i * ploidy ..
/// (i + 1) * ploidy]`.
///
/// The line comes as the source gave it, with or without its end of line,
/// which is taken off here. A line that is empty, or one whose FILTER says
/// that its variant failed a filter, gets no row at all, and the caller is
/// what leaves it out: this parses the line it is given.
///
/// # Errors
///
/// Every error of a data line of `docs/specs/io_vcf.md`, with the number of
/// the line and the column or the individual it is in. The row and the
/// genotypes are then what the parse had written when it stopped, and the
/// caller drops the block they are in.
fn parse_row(
    line: &[u8],
    number: u64,
    rules: &RowRules<'_>,
    gts: &mut [i8],
    row: &mut ParsedRow,
) -> Result<()> {
    let RowRules {
        needs,
        ploidy,
        individuals,
        ..
    } = rules;
    row.clear();
    let line = without_the_bytes_of_the_line_end(line);
    // The nine first columns are the ones whose text a block keeps, and
    // they are the ones read as text: the columns of the individuals are
    // read as bytes, which is what "Speed" of `docs/specs/io_vcf.md` asks
    // for, and no UTF-8 is checked in them.
    let (head, individual_columns) = match memchr::memchr_iter(b'\t', line).nth(8) {
        Some(ninth_tab) => (
            line.get(..ninth_tab).unwrap_or_default(),
            // The tab is at `ninth_tab`, so there is a byte after it or the
            // columns of the individuals are one empty column.
            line.get(ninth_tab.saturating_add(1)..),
        ),
        None => (line, None),
    };
    let Ok(head) = std::str::from_utf8(head) else {
        return Err(Error::VcfDataLine {
            line: number,
            place: VcfPlace::Line,
            problem: "its bytes are not valid UTF-8, and a VCF is text".to_string(),
        });
    };

    let mut columns = head.split('\t');
    let chrom_text = next_column(&mut columns, "CHROM", number)?;
    let pos_text = next_column(&mut columns, "POS", number)?;
    let id_text = next_column(&mut columns, "ID", number)?;
    let reference_text = next_column(&mut columns, "REF", number)?;
    let alternatives_text = next_column(&mut columns, "ALT", number)?;
    let quality_text = next_column(&mut columns, "QUAL", number)?;
    next_column(&mut columns, "FILTER", number)?;

    if needs.contains(Needs::CHROM_POS) {
        row.pos = parse_position(pos_text, number)?;
        row.chrom.push_str(chrom_text);
    }
    if needs.contains(Needs::ID) && id_text != MISSING_VALUE {
        row.id.push_str(id_text);
    }
    // The alleles are counted for every line that is parsed, to check the
    // allele numbers of its genotypes, also when the texts of the alleles
    // are not kept.
    row.num_alleles = count_alleles(reference_text, alternatives_text, number)?;
    if needs.contains(Needs::ALLELES) {
        row.fill_alleles(reference_text, alternatives_text);
    }
    if needs.contains(Needs::QUAL) {
        row.qual = parse_quality(quality_text, number)?.unwrap_or(f32::NAN);
    }
    // The shape of the line is checked whatever was asked for: the nine
    // first columns are there, the FORMAT has a GT key, and one column of
    // an individual comes after it at least. What is in those columns, and
    // how many of them there are, is read only when the genotypes are asked
    // for.
    //
    // INFO is not read, and its column has to be there.
    next_column(&mut columns, "INFO", number)?;
    let format_text = next_column(&mut columns, "FORMAT", number)?;
    let gt_index = gt_index_of(format_text, number)?;
    let Some(individual_columns) = individual_columns else {
        return Err(Error::VcfDataLine {
            line: number,
            place: VcfPlace::Line,
            problem: format!(
                "it has no individual column; a VCF with genotypes has the nine columns \
                 {first} and one column per individual after them",
                first = FIRST_COLUMNS.join(" "),
            ),
        });
    };
    if needs.contains(Needs::GTS) {
        fill_row_genotypes(
            gts,
            individual_columns,
            gt_index,
            individuals,
            *ploidy,
            row.num_alleles,
            number,
        )?;
    }
    Ok(())
}

/// The line without the `\n` or the `\r\n` it ends in, when it has one. The
/// genotype of the last individual is the one that would carry the `\r`.
fn without_the_bytes_of_the_line_end(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

/// The genotype of every individual of the line, into the row `gts` of the
/// block: the value of the key `GT` of each column, which an individual
/// that drops its last values still has.
///
/// The row holds one allele for each individual of the header times the
/// ploidy, and the alleles of individual `i` are `gts[i * ploidy .. (i + 1)
/// * ploidy]`.
fn fill_row_genotypes(
    gts: &mut [i8],
    individual_columns: &[u8],
    gt_index: usize,
    individuals: &[String],
    ploidy: usize,
    num_alleles: usize,
    line: u64,
) -> Result<()> {
    let expected = individuals.len().checked_mul(ploidy);
    if expected != Some(gts.len()) {
        return Err(Error::BlockArrayOfAnotherSize {
            array: "gts",
            found: gts.len(),
            expected: expected.unwrap_or(usize::MAX),
        });
    }
    let mut columns = ByteColumns::new(individual_columns);
    let mut genotypes = gts.chunks_exact_mut(ploidy);
    for (read_so_far, individual) in individuals.iter().enumerate() {
        let (Some(column), Some(genotype)) = (columns.next(), genotypes.next()) else {
            return Err(Error::VcfDataLine {
                line,
                place: VcfPlace::Line,
                problem: format!(
                    "it has the columns of {read_so_far} individuals and the header has {count}",
                    count = individuals.len(),
                ),
            });
        };
        let Some(text) = gt_of(column, gt_index) else {
            return Err(Error::VcfDataLine {
                line,
                place: VcfPlace::Individual(individual.clone()),
                problem: format!(
                    "{column} has no value where the FORMAT has GT",
                    column = shown(column),
                ),
            });
        };
        fill_row_genotype(genotype, text, num_alleles, line, individual)?;
    }
    let left_over = columns.count();
    if left_over != 0 {
        return Err(Error::VcfDataLine {
            line,
            place: VcfPlace::Line,
            problem: format!(
                "it has {left_over} {columns} more than the {count} individuals of the header",
                columns = if left_over == 1 { "column" } else { "columns" },
                count = individuals.len(),
            ),
        });
    }
    Ok(())
}

/// The bytes of the value of the key `GT` in the column of one individual,
/// which is the `gt_index`th of the values the column holds, and `None`
/// when the column has fewer values than that.
fn gt_of(column: &[u8], gt_index: usize) -> Option<&[u8]> {
    let mut value = column;
    for _ in 0..gt_index {
        let at = memchr::memchr(b':', value)?;
        value = value.get(at.saturating_add(1)..)?;
    }
    match memchr::memchr(b':', value) {
        Some(at) => value.get(..at),
        None => Some(value),
    }
}

/// The alleles of the genotype of one individual, into the `ploidy` alleles
/// of that individual in the row of the block.
///
/// A genotype written as a single dot is a missing genotype of the ploidy of
/// the file, and any other number of alleles than the ploidy is an error,
/// which the length of `genotype` is what says: it holds the ploidy the
/// reader was asked for.
fn fill_row_genotype(
    genotype: &mut [i8],
    text: &[u8],
    num_alleles: usize,
    line: u64,
    individual: &str,
) -> Result<()> {
    // VCF 4.4 lets a genotype start with its separator, `/0/1`.
    let text = match text.first() {
        Some(b'/' | b'|') => text.get(1..).unwrap_or_default(),
        _ => text,
    };
    if text == MISSING_VALUE.as_bytes() {
        for allele in genotype.iter_mut() {
            *allele = MISSING_ALLELE;
        }
        return Ok(());
    }
    let mut written: usize = 0;
    let mut alleles = genotype.iter_mut();
    for allele_text in text.split(|byte| *byte == b'/' || *byte == b'|') {
        let allele = parse_row_allele(allele_text, num_alleles, line, individual)?;
        if let Some(place) = alleles.next() {
            *place = allele;
        }
        // The alleles of a genotype are at most the bytes of its text, so
        // the count does not reach the largest `usize`.
        written = written.saturating_add(1);
    }
    if written != genotype.len() {
        return Err(Error::VcfGenotypePloidy {
            line,
            individual: individual.to_string(),
            found: written,
            expected: genotype.len(),
        });
    }
    Ok(())
}

/// One allele of a genotype, as bytes: a number of the alleles the variant
/// declares, or [`MISSING_ALLELE`] for a dot.
///
/// The bytes are never turned into text, so a byte that is not text is a
/// byte that is not a digit, and the message shows it as the replacement
/// character.
fn parse_row_allele(text: &[u8], num_alleles: usize, line: u64, individual: &str) -> Result<i8> {
    if text == MISSING_VALUE.as_bytes() {
        return Ok(MISSING_ALLELE);
    }
    let wrong = |problem: String| Error::VcfDataLine {
        line,
        place: VcfPlace::Individual(individual.to_string()),
        problem,
    };
    if text.is_empty() {
        return Err(wrong("`` is not an allele number".to_string()));
    }
    let mut number: u32 = 0;
    for byte in text {
        let Some(digit) = char::from(*byte).to_digit(10) else {
            return Err(wrong(format!(
                "{text} is not an allele number, which is a run of digits",
                text = shown(text),
            )));
        };
        // Once the number is above the largest allele the answer is the
        // same whatever its other digits are, so it stops growing there and
        // the two operations cannot overflow.
        number = number.saturating_mul(10).saturating_add(digit);
    }
    let Ok(allele) = i8::try_from(number) else {
        return Err(wrong(format!(
            "the allele {text} is above {MAX_ALLELE}, the largest allele popnei holds",
            text = shown(text),
        )));
    };
    // The allele is not negative, since its text was parsed as a `u32`.
    if usize::from(allele.unsigned_abs()) >= num_alleles {
        return Err(wrong(format!(
            "the allele {number} is not one of the {num_alleles} alleles that REF and ALT declare"
        )));
    }
    Ok(allele)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::{BufRead, BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{
        BYTES_PER_BATCH, BatchRow, GZIP_FLAGS, LINES_PER_BATCH, MAX_PLOIDY, MISSING_VALUE,
        ParsedRow, RowRules, VcfOptions, VcfPlace, VcfReader, WithTheLastBytes, parse_row,
        parse_rows, parse_rows_one_by_one, written_by_bgzip,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::{Error, Result};
    use crate::variant::{MISSING_ALLELE, Needs};

    /// The panic that the test of a reader whose parse did not come back
    /// injects into the parse of one line. It is the only way to make the
    /// parse panic: no VCF does it, and every error of a line is a value
    /// that the line carries back.
    pub(super) fn panic_if_the_test_asked_for_it(number: u64, rules: &super::RowRules<'_>) {
        assert!(
            rules.panic_at_line != Some(number),
            "the parse of the line {number} panicked, which this test asked for"
        );
    }

    /// The reference VCFs and what bcftools 1.24 read in them live at the
    /// root of the repository, beside the Python tests that read the same
    /// files, and not inside this crate. The path is built from the
    /// directory of the manifest, so it holds whether the tests are run
    /// with `cargo test --workspace` or with `cargo test -p popnei`.
    fn reference_vcf(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/vcf")
            .join(name)
    }

    /// The error of a VCF that `VcfReader::new` has to refuse.
    fn error_of(vcf: &str, options: VcfOptions) -> Error {
        match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => panic!(
                "the reader was built, over the individuals {:?}",
                reader.individuals()
            ),
            Err(error) => error,
        }
    }

    /// The header of the four files of the reference, with three
    /// individuals.
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

    /// The first data line of a VCF that `vcf_of` writes is the fourth line
    /// of the file: the errors count the three lines of the header.
    const FIRST_DATA_LINE: u64 = 4;

    /// The options of a file of another ploidy, or of one read with every
    /// variant given, with the size of a block that popnei chooses.
    fn options(ploidy: usize, only_passed: bool) -> VcfOptions {
        VcfOptions {
            ploidy,
            only_passed,
            num_vars_per_block: None,
        }
    }

    /// The same options with the blocks of the size that a test asks for,
    /// which is how a test reads a file of a few lines in several blocks.
    fn in_blocks_of(options: VcfOptions, num_vars_per_block: usize) -> VcfOptions {
        VcfOptions {
            num_vars_per_block: Some(num_vars_per_block),
            ..options
        }
    }

    /// A reader over bytes held in memory, which is how the tests give a
    /// VCF of their own.
    fn reader_over(vcf: &str, options: VcfOptions) -> VcfReader<Cursor<Vec<u8>>> {
        match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

    /// A reader over one of the reference files.
    fn reader_of_file(name: &str, options: VcfOptions) -> VcfReader<BufReader<File>> {
        match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        }
    }

    /// One variant as the tables of "How it is verified" of
    /// `docs/specs/io_vcf.md` give it, with the name of the chromosome in
    /// the place of the number the reader gave it. A field the block does
    /// not hold is empty, 0 or `None`.
    #[derive(Debug, Clone, PartialEq)]
    struct Row {
        chrom: String,
        pos: u64,
        id: String,
        alleles: Vec<String>,
        qual: Option<f32>,
        gts: Vec<i8>,
    }

    fn row(
        chrom: &str,
        pos: u64,
        id: &str,
        alleles: &[&str],
        qual: Option<f32>,
        gts: &[i8],
    ) -> Row {
        Row {
            chrom: chrom.to_string(),
            pos,
            id: id.to_string(),
            alleles: alleles.iter().map(|text| (*text).to_string()).collect(),
            qual,
            gts: gts.to_vec(),
        }
    }

    /// Every block a reader gives, until it has no more or it fails. Each
    /// of them holds one variant at least and is of its size, which the
    /// trait of a reader of blocks asks of every reader.
    fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
            assert!(block.num_vars > 0, "a block of no variant");
            block.check().expect("the block is of its size");
            assert_eq!(block.num_individuals, reader.individuals().len());
            assert_eq!(block.ploidy, reader.ploidy());
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// The variants of the blocks of a reader, joined, through the views of
    /// one variant, with the names of their chromosomes.
    fn rows_of(reader: &mut impl BlockReader) -> Result<Vec<Row>> {
        let mut rows = Vec::new();
        loop {
            let Some(block) = reader.next_block()? else {
                return Ok(rows);
            };
            assert!(block.num_vars > 0, "a block of no variant");
            block.check().expect("the block is of its size");
            for view in block.variants() {
                rows.push(Row {
                    chrom: view
                        .chrom()
                        .and_then(|number| reader.chroms().name(number))
                        .unwrap_or_default()
                        .to_string(),
                    pos: view.pos().unwrap_or_default(),
                    id: view.id().unwrap_or_default().to_string(),
                    alleles: (0..view.num_alleles().unwrap_or_default())
                        .map(|allele| view.allele(allele).unwrap_or_default().to_string())
                        .collect(),
                    qual: view.qual().filter(|qual| !qual.is_nan()),
                    gts: view.gts().to_vec(),
                });
            }
        }
    }

    /// The variants of a VCF written in a test.
    fn rows_read(vcf: &str, options: VcfOptions) -> Vec<Row> {
        match rows_of(&mut reader_over(vcf, options)) {
            Ok(rows) => rows,
            Err(error) => panic!("the reader stopped at {error}"),
        }
    }

    /// The variants of one of the reference files.
    fn rows_of_file(name: &str, options: VcfOptions) -> Vec<Row> {
        match rows_of(&mut reader_of_file(name, options)) {
            Ok(rows) => rows,
            Err(error) => panic!("{name}: the reader stopped at {error}"),
        }
    }

    /// The error a VCF written in a test stops the reader at.
    fn error_reading(vcf: &str, options: VcfOptions) -> Error {
        match rows_of(&mut reader_over(vcf, options)) {
            Ok(rows) => panic!("the reader gave {} variants and no error", rows.len()),
            Err(error) => error,
        }
    }

    /// The four rows of the table of `cases.vcf` of "How it is verified" of
    /// `docs/specs/io_vcf.md`, which bcftools 1.24 printed.
    fn the_rows_of_cases() -> Vec<Row> {
        vec![
            row(
                "chr1",
                100,
                "rs1",
                &["A", "T"],
                Some(29.5),
                &[0, 0, 0, 1, 1, 1],
            ),
            row(
                "chr1",
                200,
                "",
                &["A", "T"],
                None,
                &[MISSING_ALLELE, MISSING_ALLELE, 0, 1, MISSING_ALLELE, 0],
            ),
            row(
                "chr1",
                300,
                "",
                &["A", "G", "T"],
                Some(67.0),
                &[1, 2, 2, 1, 2, 2],
            ),
            row("chr1", 400, "", &["T"], Some(47.0), &[0, 0, 0, 0, 0, 0]),
        ]
    }

    /// The two rows of the table of `differences.vcf` of the same section,
    /// the two variants that pyNei does not read as bcftools does.
    fn the_rows_of_differences() -> Vec<Row> {
        vec![
            row(
                "chr2",
                50,
                "ms1",
                &["GTC", "G", "GTCT"],
                Some(50.0),
                &[0, 1, 0, 2, MISSING_ALLELE, MISSING_ALLELE],
            ),
            row(
                "chr2",
                60,
                "",
                &["A", "<DEL>", "*"],
                None,
                &[0, 1, 2, 2, 0, 0],
            ),
        ]
    }

    /// The batches to read a file of a few lines in, so that a test of what
    /// the reader does at a line is made with that line in a batch of its
    /// own, at the start of a batch, in the middle of one and with the
    /// whole file in one: the batching is what carries the lines and the
    /// errors of a file to the blocks that give them.
    const BATCHES_TO_TRY: [usize; 4] = [LINES_PER_BATCH, 4, 2, 1];

    // What the reader reads in a header, and what it refuses there.

    #[test]
    fn the_individuals_of_a_plain_and_of_a_gzipped_vcf_are_read() {
        for name in ["cases.vcf", "cases.vcf.gz"] {
            let reader = reader_of_file(name, VcfOptions::default());
            assert_eq!(reader.individuals(), ["ind1", "ind2", "ind3"], "{name}");
            assert_eq!(reader.ploidy(), 2, "{name}");
            assert!(reader.chroms().is_empty(), "{name}");
        }
    }

    #[test]
    fn the_gzip_bytes_are_found_in_a_source_that_is_not_a_file() {
        let bytes = std::fs::read(reference_vcf("differences.vcf.gz")).unwrap();
        assert_eq!(bytes.get(..2), Some([0x1f, 0x8b].as_slice()));
        let reader = VcfReader::new(Cursor::new(bytes), VcfOptions::default()).unwrap();
        assert_eq!(reader.individuals(), ["ind1", "ind2", "ind3"]);
    }

    /// A source whose buffer holds one byte at a time, which is what a pipe
    /// can give and what a page that hands over the bytes of a file a few
    /// at a time would be: the two bytes of gzip are never in its buffer
    /// together.
    fn one_byte_at_a_time(bytes: Vec<u8>) -> BufReader<Cursor<Vec<u8>>> {
        BufReader::with_capacity(1, Cursor::new(bytes))
    }

    #[test]
    fn a_path_that_no_file_is_at_gives_an_error_that_carries_the_path() {
        let path = reference_vcf("no_such_file.vcf");
        let error = match VcfReader::from_path(&path, VcfOptions::default()) {
            Ok(reader) => panic!("the reader was built over {:?}", reader.individuals()),
            Err(error) => error,
        };
        let Error::FileNotOpened {
            path: named,
            source,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(named, path);
        assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
        let message = Error::FileNotOpened {
            path: named,
            source,
        }
        .to_string();
        assert!(message.contains("no_such_file.vcf"), "{message}");
    }

    #[test]
    fn a_source_that_gives_one_byte_at_a_time_is_read_gzipped_and_plain() {
        for name in ["cases.vcf", "cases.vcf.gz"] {
            let bytes = std::fs::read(reference_vcf(name)).unwrap();
            let mut reader = VcfReader::new(one_byte_at_a_time(bytes), options(2, false)).unwrap();
            assert_eq!(reader.individuals(), ["ind1", "ind2", "ind3"], "{name}");
            assert_eq!(rows_of(&mut reader).unwrap(), the_rows_of_cases(), "{name}");
        }
    }

    #[test]
    fn a_source_that_starts_with_neither_a_hash_nor_the_gzip_bytes_is_refused() {
        let error = error_of("chr1\t100\trs1\n", VcfOptions::default());
        let Error::NotAVcf { found } = error else {
            panic!("the error is {error}");
        };
        assert!(found.contains("chr1"), "{found}");
    }

    #[test]
    fn a_header_with_no_format_column_is_refused() {
        let vcf = "##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tind1\tind2
";
        let error = error_of(vcf, VcfOptions::default());
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("FORMAT"), "{problem}");
    }

    #[test]
    fn a_header_with_a_format_column_and_no_individual_is_refused() {
        let vcf = "##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT
";
        let error = error_of(vcf, VcfOptions::default());
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("individual"), "{problem}");
    }

    #[test]
    fn two_individuals_with_the_same_name_are_refused() {
        let vcf = "##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind1
";
        let error = error_of(vcf, VcfOptions::default());
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("ind1"), "{problem}");
    }

    #[test]
    fn an_individual_with_no_name_is_refused() {
        let vcf = "##fileformat=VCFv4.4
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\t
";
        let error = error_of(vcf, VcfOptions::default());
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("no name"), "{problem}");
    }

    #[test]
    fn a_header_that_ends_before_the_chrom_line_is_refused() {
        let vcf = "##fileformat=VCFv4.4\n##contig=<ID=chr1>\n";
        let error = error_of(vcf, VcfOptions::default());
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        // The other errors of a header name the #CHROM line too, so what
        // this one has to say is that there is none.
        assert!(problem.contains("no #CHROM line"), "{problem}");
    }

    #[test]
    fn a_header_line_whose_bytes_are_not_text_is_refused() {
        let mut vcf = b"##fileformat=VCFv4.4\n##contig=<ID=\xffchr1>\n".to_vec();
        vcf.extend_from_slice(HEADER.as_bytes());
        let error = match VcfReader::new(Cursor::new(vcf), VcfOptions::default()) {
            Ok(reader) => panic!("the reader was built over {:?}", reader.individuals()),
            Err(error) => error,
        };
        let Error::VcfHeader { problem } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("UTF-8"), "{problem}");
        assert!(problem.contains('2'), "{problem}");
    }

    #[test]
    fn an_empty_source_is_refused() {
        let error = error_of("", VcfOptions::default());
        let Error::NotAVcf { found } = error else {
            panic!("the error is {error}");
        };
        assert!(found.contains("no byte"), "{found}");
    }

    #[test]
    fn a_gzipped_source_that_is_not_a_vcf_is_refused() {
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, b"hello\n").unwrap();
        let gzipped = encoder.finish().unwrap();
        let error = match VcfReader::new(Cursor::new(gzipped), VcfOptions::default()) {
            Ok(reader) => panic!("the reader was built over {:?}", reader.individuals()),
            Err(error) => error,
        };
        let Error::NotAVcf { found } = error else {
            panic!("the error is {error}");
        };
        assert!(found.contains("hello"), "{found}");
    }

    #[test]
    fn a_ploidy_of_zero_and_one_above_the_largest_are_refused() {
        for asked_for in [0, MAX_PLOIDY.saturating_add(1), usize::MAX] {
            let error = error_of(HEADER, options(asked_for, true));
            let Error::VcfPloidyOutOfRange { ploidy, largest } = error else {
                panic!("the error of the ploidy {asked_for} is {error}");
            };
            assert_eq!((ploidy, largest), (asked_for, MAX_PLOIDY));
        }
    }

    #[test]
    fn a_reader_of_blocks_of_no_variant_is_refused() {
        // A block holds one variant at least, and the caller that wants the
        // size popnei chooses asks for none instead of asking for 0.
        let error = error_of(HEADER, in_blocks_of(VcfOptions::default(), 0));
        assert!(matches!(error, Error::BlockOfNoVariants), "{error}");
        let message = error.to_string();
        assert!(message.contains('0'), "{message}");
    }

    /// The genotypes of a block are its variants times the individuals
    /// times the ploidy, and a size that a caller wrote carries that
    /// multiplication beyond what a `usize` holds. It is an error when the
    /// reader is built, before a line of the file is read, and not a panic,
    /// on a machine of 32 bit addresses as on one of 64.
    #[test]
    fn a_block_of_more_genotypes_than_a_usize_holds_is_refused_when_the_reader_is_built() {
        let error = error_of(HEADER, in_blocks_of(VcfOptions::default(), usize::MAX));
        let message = error.to_string();
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
        // The size is the one the caller wrote, so what the message says to
        // do about it is to write a smaller one.
        assert!(message.contains("ask for fewer variants"), "{message}");
        assert!(!message.contains("popnei chose"), "{message}");
    }

    /// The size that popnei chooses is not checked when the reader is
    /// built: a caller that opens a file for its individuals alone, which
    /// is what `open_vcf` and `openVcf` do, is not refused for a size that
    /// nobody asked for. A `usize` is 64 bits natively, and the block of
    /// 100 variants of 170000 individuals of the ploidy 255 that is more
    /// than one holds in wasm is 4.3 thousand million here, so what this
    /// test sees is that the reader is built and gives the variants; the
    /// node test of `js/popnei` is where the refusal would show.
    #[test]
    fn the_size_popnei_chooses_is_not_refused_when_the_reader_is_built() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, options(MAX_PLOIDY, true));
        assert_eq!(reader.individuals().len(), 3);
        // The genotypes of that line are not of the ploidy 255, so what is
        // asserted is the error of the line and not one of the size.
        let error = match rows_of(&mut reader) {
            Ok(rows) => panic!("the reader gave {} variants", rows.len()),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::VcfGenotypePloidy { .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn the_largest_ploidy_is_read() {
        let mut genotype = String::from("0");
        for _ in 1..MAX_PLOIDY {
            genotype.push_str("/1");
        }
        let line = format!("chr1 100 . A T . PASS . GT {genotype} . .");
        let rows = rows_read(&vcf_of(&[&line]), options(MAX_PLOIDY, true));
        let first = rows.first().expect("one variant");
        assert_eq!(first.gts.len(), MAX_PLOIDY.saturating_mul(3));
        assert_eq!(first.gts.first(), Some(&0));
        assert_eq!(first.gts.last(), Some(&MISSING_ALLELE));
    }

    // What a data line gives, and what is refused in one. The tests of the
    // row parser above are made at one line; these are made at the blocks
    // the reader gives.

    #[test]
    fn the_four_variants_of_cases_vcf_are_read_when_every_variant_is_given() {
        for name in ["cases.vcf", "cases.vcf.gz"] {
            assert_eq!(
                rows_of_file(name, options(2, false)),
                the_rows_of_cases(),
                "{name}"
            );
        }
    }

    #[test]
    fn the_default_leaves_out_the_variant_of_cases_vcf_that_failed_its_filter() {
        let all = the_rows_of_cases();
        // The second variant, chr1 200, is the one whose FILTER is q10.
        let passed = vec![all[0].clone(), all[2].clone(), all[3].clone()];
        for name in ["cases.vcf", "cases.vcf.gz"] {
            assert_eq!(rows_of_file(name, VcfOptions::default()), passed, "{name}");
        }
    }

    #[test]
    fn the_leading_separators_and_the_dot_of_differences_vcf_are_read() {
        for only_passed in [true, false] {
            for name in ["differences.vcf", "differences.vcf.gz"] {
                assert_eq!(
                    rows_of_file(name, options(2, only_passed)),
                    the_rows_of_differences(),
                    "{name}, only_passed {only_passed}"
                );
            }
        }
    }

    #[test]
    fn the_chromosomes_are_numbered_in_the_order_of_the_variants_that_are_given() {
        let vcf = vcf_of(&[
            "chr9 10 . A T . q10 . GT 0/0 0/1 1/1",
            "chr1 20 rs2 A T 9.5 PASS . GT 0/0 0/1 1/1",
            "chr9 30 . A T . PASS . GT 0/0 0/1 1/1",
        ]);
        // The variant that the FILTER leaves out is the first of the file
        // and its chromosome gets no number, so `chr1` is the number 0
        // although `chr9` comes first in the file. The blocks of one
        // variant are what makes the numbers of a file read in several
        // blocks the same as in one.
        for num_vars_per_block in [1, 2, 100] {
            let mut reader = reader_over(
                &vcf,
                in_blocks_of(VcfOptions::default(), num_vars_per_block),
            );
            let rows = rows_of(&mut reader).unwrap();
            assert_eq!(
                rows.iter().map(|row| row.pos).collect::<Vec<u64>>(),
                [20, 30],
                "blocks of {num_vars_per_block}"
            );
            assert_eq!(
                rows.iter()
                    .map(|row| row.chrom.as_str())
                    .collect::<Vec<&str>>(),
                ["chr1", "chr9"],
                "blocks of {num_vars_per_block}"
            );
            assert_eq!(reader.chroms().len(), 2);
            assert_eq!(reader.chroms().name(0), Some("chr1"));
            assert_eq!(reader.chroms().name(1), Some("chr9"));
        }
    }

    #[test]
    fn a_tetraploid_genotype_with_a_ploidy_of_two_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0/1/1 0/1/1/1 0/0/0/0"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfGenotypePloidy {
            line,
            individual,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (line, individual.as_str(), found, expected),
            (FIRST_DATA_LINE, "ind1", 4, 2)
        );
    }

    #[test]
    fn a_haploid_genotype_with_a_ploidy_of_two_is_refused_after_the_blocks_before_it() {
        let vcf = vcf_of(&[
            "chr1 100 . A T . PASS . GT 0/0 0/1 1/1",
            "chr1 200 . A T . PASS . GT 1 0 1",
        ]);
        let mut reader = reader_over(&vcf, in_blocks_of(VcfOptions::default(), 1));

        let first = reader.next_block().unwrap().expect("the first block");
        assert_eq!(
            (first.num_vars, first.pos.as_deref()),
            (1, Some([100].as_slice()))
        );

        let error = reader.next_block().unwrap_err();
        let Error::VcfGenotypePloidy {
            line,
            individual,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (line, individual.as_str(), found, expected),
            (5, "ind1", 1, 2)
        );
        // The block the wrong line was in is lost and there is no block
        // after it.
        assert!(reader.next_block().unwrap().is_none());
    }

    #[test]
    fn a_tetraploid_vcf_read_with_a_ploidy_of_four_is_read() {
        let vcf = vcf_of(&[
            "chr1 100 . A T . PASS . GT 0/0/1/1 0/1/1/1 0/0/0/0",
            "chr1 200 . A T . PASS . GT . . .",
        ]);
        let missing = [MISSING_ALLELE; 12];
        assert_eq!(
            rows_read(&vcf, options(4, true)),
            vec![
                row(
                    "chr1",
                    100,
                    "",
                    &["A", "T"],
                    None,
                    &[0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0]
                ),
                row("chr1", 200, "", &["A", "T"], None, &missing),
            ]
        );
    }

    #[test]
    fn an_allele_that_the_variant_does_not_declare_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/2 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        // The alleles of ALT are counted for every line that is parsed, so
        // that the allele numbers are checked also when the texts of the
        // alleles are not kept.
        reader.set_needs(Needs::GTS);

        let error = reader.next_block().unwrap_err();
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(line, FIRST_DATA_LINE);
        assert_eq!(place, VcfPlace::Individual("ind2".to_string()));
        assert!(problem.contains('2'), "{problem}");
        assert!(problem.contains("declare"), "{problem}");
    }

    #[test]
    fn an_allele_above_the_largest_one_popnei_holds_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/128 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(line, FIRST_DATA_LINE);
        assert_eq!(place, VcfPlace::Individual("ind2".to_string()));
        assert!(problem.contains("128"), "{problem}");
        assert!(problem.contains("127"), "{problem}");
    }

    #[test]
    fn a_line_with_the_genotypes_of_two_individuals_under_a_header_of_three_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
        assert!(problem.contains('2') && problem.contains('3'), "{problem}");
    }

    #[test]
    fn a_line_with_the_genotypes_of_four_individuals_under_a_header_of_three_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1 0/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
        assert!(problem.contains('3'), "{problem}");
        assert!(problem.contains("1 column more"), "{problem}");

        // Two columns more, where the count is plural.
        let vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1 0/1 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine { problem, .. } = error else {
            panic!("the error is {error}");
        };
        assert!(problem.contains("2 columns more"), "{problem}");
    }

    #[test]
    fn a_format_with_no_gt_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . DP 3 4 5"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("FORMAT")));
        assert!(problem.contains("GT"), "{problem}");
    }

    #[test]
    fn a_position_that_is_not_a_number_is_refused() {
        let vcf = vcf_of(&["chr1 x . A T . PASS . GT 0/0 0/1 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("POS")));
        assert!(problem.contains('x'), "{problem}");
    }

    /// A column that is not parsed is not checked, which "How it runs" of
    /// `docs/specs/io_vcf.md` decides: the reader that this one took the
    /// place of parsed the position of every line, and a position that is
    /// not a number was an error whatever was asked for.
    #[test]
    fn a_position_that_is_not_a_number_is_read_when_the_position_was_not_asked_for() {
        let vcf = vcf_of(&["chr1 x . A T . PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::GTS);
        let rows = rows_of(&mut reader).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].gts, [0, 0, 0, 1, 1, 1]);
        assert_eq!(rows[0].pos, 0);
    }

    #[test]
    fn a_quality_that_is_not_a_number_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T x PASS . GT 0/0 0/1 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("QUAL")));
        assert!(problem.contains('x'), "{problem}");
    }

    #[test]
    fn an_allele_number_is_a_run_of_digits_and_nothing_else() {
        // `+1` and `-1` are numbers that Rust's own parser reads, and the
        // VCF has neither; an allele of no digit at all is not one either.
        for genotype in ["0/+1", "0/-1", "0/1x", "0/", "/", "0/1.5"] {
            let line = format!("chr1 100 . A T . PASS . GT 0/0 {genotype} 1/1");
            let error = error_reading(&vcf_of(&[&line]), VcfOptions::default());
            let Error::VcfDataLine { line, place, .. } = error else {
                panic!("the error of `{genotype}` is {error}");
            };
            assert_eq!(
                (line, place),
                (FIRST_DATA_LINE, VcfPlace::Individual("ind2".to_string())),
                "{genotype}"
            );
        }
    }

    #[test]
    fn an_allele_number_that_no_i8_holds_is_refused_for_being_above_the_largest() {
        // 4294967296 does not fit in the u32 the number was parsed into,
        // and the answer is the same as for 128: it is above 127.
        for genotype in ["0/128", "0/4294967296", "0/99999999999999999999"] {
            let line = format!("chr1 100 . A T . PASS . GT 0/0 {genotype} 1/1");
            let error = error_reading(&vcf_of(&[&line]), VcfOptions::default());
            let Error::VcfDataLine { problem, .. } = error else {
                panic!("the error of `{genotype}` is {error}");
            };
            assert!(problem.contains("127"), "{genotype}: {problem}");
        }
    }

    #[test]
    fn a_gt_that_is_not_the_first_key_of_the_format_is_read() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . DP:GT 3:0/1 4:1/1 5:./."]);
        assert_eq!(
            rows_read(&vcf, VcfOptions::default()),
            vec![row(
                "chr1",
                100,
                "",
                &["A", "T"],
                None,
                &[0, 1, 1, 1, MISSING_ALLELE, MISSING_ALLELE]
            )]
        );
    }

    #[test]
    fn a_vcf_whose_lines_end_in_a_carriage_return_is_read() {
        let vcf = "##fileformat=VCFv4.4\r\n\
                   #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3\r\n\
                   chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t.|.\r\n";
        assert_eq!(
            rows_read(vcf, VcfOptions::default()),
            vec![row(
                "chr1",
                100,
                "",
                &["A", "T"],
                None,
                &[0, 0, 0, 1, MISSING_ALLELE, MISSING_ALLELE]
            )]
        );
    }

    #[test]
    fn an_empty_line_at_the_end_is_skipped() {
        let vcf = format!("{}\n", vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]));
        assert_eq!(
            rows_read(&vcf, VcfOptions::default()),
            vec![row("chr1", 100, "", &["A", "T"], None, &[0, 0, 0, 1, 1, 1])]
        );
    }

    #[test]
    fn a_last_line_with_no_end_of_line_is_read() {
        let vcf = format!(
            "{HEADER}chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
             chr1\t200\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
        );
        assert_eq!(
            rows_read(&vcf, VcfOptions::default()),
            vec![
                row("chr1", 100, "", &["A", "T"], None, &[0, 0, 0, 1, 1, 1]),
                row("chr1", 200, "", &["A", "T"], None, &[0, 0, 0, 1, 1, 1]),
            ]
        );
    }

    #[test]
    fn a_vcf_with_a_header_and_no_variant_gives_no_block_and_no_error() {
        let mut reader = reader_over(HEADER, VcfOptions::default());
        assert!(reader.next_block().unwrap().is_none());
        assert!(reader.next_block().unwrap().is_none());
        assert!(reader.chroms().is_empty());
    }

    #[test]
    fn the_blocks_after_the_error_of_the_first_line_of_a_file_are_no_blocks() {
        // The wrong line first and three good ones after it, in blocks of
        // one variant: a reader that did not end at its error would give
        // the three, which is what `docs/specs/block.md` refuses of every
        // reader.
        let vcf = vcf_of(&[
            "chr1 100 . A T . PASS . GT 0/0 0/1 1",
            "chr1 200 . A T . PASS . GT 0/0 0/1 1/1",
            "chr1 300 . A T . PASS . GT 0/0 0/1 1/1",
            "chr1 400 . A T . PASS . GT 0/0 0/1 1/1",
        ]);
        let mut reader = reader_over(&vcf, in_blocks_of(VcfOptions::default(), 1));

        let error = match reader.next_block() {
            Ok(block) => panic!("the reader gave {block:?}"),
            Err(error) => error,
        };
        assert!(
            matches!(error, Error::VcfGenotypePloidy { line: 4, .. }),
            "the error is {error}"
        );
        for call in 1..=3 {
            assert!(
                reader.next_block().expect("no block").is_none(),
                "the call {call} after the error"
            );
        }
    }

    #[test]
    fn a_block_after_an_error_and_after_the_last_block_is_no_block() {
        let mut reader = reader_over(
            &vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]),
            VcfOptions::default(),
        );
        assert!(reader.next_block().unwrap().is_some());
        assert!(reader.next_block().unwrap().is_none());
        assert!(reader.next_block().unwrap().is_none());

        let mut reader = reader_over(
            &vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1"]),
            VcfOptions::default(),
        );
        assert!(reader.next_block().is_err());
        assert!(reader.next_block().unwrap().is_none());
        assert!(reader.next_block().unwrap().is_none());
    }

    #[test]
    fn an_allele_that_alt_declares_and_no_genotype_carries_is_read() {
        let vcf = vcf_of(&["chr1 100 . A G,T . PASS . GT 0/1 0/0 1/1"]);
        assert_eq!(
            rows_read(&vcf, VcfOptions::default()),
            vec![row(
                "chr1",
                100,
                "",
                &["A", "G", "T"],
                None,
                &[0, 1, 0, 0, 1, 1]
            )]
        );
    }

    #[test]
    fn a_line_that_failed_its_filter_is_skipped_before_its_genotypes_are_read() {
        let vcf = vcf_of(&["chr1 100 . A T . q10 . GT 0/0/1/1 0/1 0/1"]);
        assert!(rows_read(&vcf, VcfOptions::default()).is_empty());

        let error = error_reading(&vcf, options(2, false));
        assert!(
            matches!(error, Error::VcfGenotypePloidy { found: 4, .. }),
            "the error is {error}"
        );
    }

    #[test]
    fn a_line_that_ends_before_its_filter_is_refused_with_the_default() {
        // Whether this line would be skipped cannot be known: the FILTER
        // is the seventh column and the line has four. The columns up to
        // the FILTER have to be there, and only what is inside them is
        // read late.
        let vcf = vcf_of(&["chr1 10 . A"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
    }

    #[test]
    fn a_line_that_failed_its_filter_is_not_read_before_its_filter_either() {
        let vcf = vcf_of(&[
            "chr9 x . A T . q10 . GT 0/0 0/1 1/1",
            "chr1 200 . A T . PASS . GT 0/0 0/1 1/1",
        ]);
        assert_eq!(
            rows_read(&vcf, VcfOptions::default()),
            vec![row("chr1", 200, "", &["A", "T"], None, &[0, 0, 0, 1, 1, 1])]
        );
    }

    #[test]
    fn with_the_genotypes_alone_a_block_has_no_column() {
        let vcf = vcf_of(&["chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::GTS);

        let block = reader.next_block().unwrap().expect("a block");
        assert_eq!(block.num_vars, 1);
        assert_eq!(block.gts, [0, 0, 0, 1, 1, 1]);
        assert!(block.chrom.is_none());
        assert!(block.pos.is_none());
        assert!(block.id.is_none());
        assert!(block.alleles.is_none());
        assert!(block.qual.is_none());
        assert_eq!(block.fields(), Needs::GTS);
        // The chromosome of a line that was not asked for its chromosome
        // gets no number.
        assert!(reader.chroms().is_empty());
    }

    #[test]
    fn with_the_id_and_the_alleles_alone_a_block_has_those_two_columns_and_no_genotype() {
        let vcf = vcf_of(&["chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::ID | Needs::ALLELES);

        let block = reader.next_block().unwrap().expect("a block");
        assert_eq!(block.num_vars, 1);
        assert!(block.gts.is_empty());
        assert_eq!(block.id.as_deref(), Some(["rs1".to_string()].as_slice()));
        let alleles = block.alleles.as_ref().expect("the alleles");
        assert_eq!(alleles.num_alleles(0), 2);
        assert_eq!(alleles.allele(0, 1), "T");
        assert!(block.chrom.is_none());
        assert!(block.qual.is_none());
        assert_eq!(block.fields(), Needs::ID | Needs::ALLELES);
    }

    /// The variants of a VCF read with the id and the alleles asked for
    /// and not the genotypes, or the error the reader stops at.
    fn read_without_the_genotypes(vcf: &str) -> Result<Vec<Row>> {
        let mut reader = reader_over(vcf, VcfOptions::default());
        reader.set_needs(Needs::ID | Needs::ALLELES);
        rows_of(&mut reader)
    }

    #[test]
    fn a_line_of_seven_columns_is_refused_with_the_genotypes_not_asked_for() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS"]);
        let error = read_without_the_genotypes(&vcf).unwrap_err();
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
    }

    #[test]
    fn a_format_with_no_gt_is_refused_with_the_genotypes_not_asked_for() {
        let vcf = vcf_of(&["chr1 100 . A T . PASS . DP 3 4 5"]);
        let error = read_without_the_genotypes(&vcf).unwrap_err();
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("FORMAT")));
    }

    #[test]
    fn what_is_in_the_columns_of_the_individuals_is_not_read_without_the_genotypes() {
        // A tetraploid genotype under a ploidy of 2, and a line with the
        // columns of two individuals under a header with three: errors
        // when the genotypes are asked for, and not looked at here.
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0/1/1 0/1 0/1",
            "chr1 200 rs2 A T . PASS . GT 0/0 0/1",
        ]);
        assert_eq!(
            read_without_the_genotypes(&vcf).unwrap(),
            // The position was not asked for either, so it is not parsed
            // and the column of a block does not hold it.
            vec![
                row("", 0, "rs1", &["A", "T"], None, &[]),
                row("", 0, "rs2", &["A", "T"], None, &[]),
            ]
        );
    }

    #[test]
    fn a_reader_asked_for_the_genotypes_alone_gives_its_next_block_without_the_columns() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T 29.5 PASS . GT 0/0 0/1 1/1",
        ]);
        let mut reader = reader_over(&vcf, in_blocks_of(VcfOptions::default(), 1));

        let first = reader.next_block().unwrap().expect("the first block");
        assert_eq!(first.fields(), Needs::ALL);
        assert_eq!(first.id.as_deref(), Some(["rs1".to_string()].as_slice()));

        reader.set_needs(Needs::GTS);
        let second = reader.next_block().unwrap().expect("the second block");
        assert_eq!(second.fields(), Needs::GTS);
        assert!(second.id.is_none());
        assert!(second.alleles.is_none());
        assert!(second.qual.is_none());
        assert_eq!(second.gts, [0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn a_data_line_whose_bytes_are_not_text_is_refused_with_its_number() {
        let mut vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]).into_bytes();
        vcf.extend_from_slice(b"chr1\t200\t.\t\xffA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n");
        for lines_per_batch in BATCHES_TO_TRY {
            let mut reader = VcfReader::new(
                Cursor::new(vcf.clone()),
                in_blocks_of(VcfOptions::default(), 1),
            )
            .unwrap();
            reader.set_lines_per_batch(lines_per_batch);

            // The line that cannot be read is the one after a good one,
            // whatever batch they fall in: the block of the good line is
            // given first and the error comes at the call after it.
            let first = reader
                .next_block()
                .unwrap_or_else(|error| panic!("{lines_per_batch}: {error}"))
                .expect("the first block");
            assert_eq!(first.pos.as_deref(), Some([100].as_slice()));

            let error = reader.next_block().unwrap_err();
            let Error::VcfDataLine {
                line,
                place,
                problem,
            } = error
            else {
                panic!("in batches of {lines_per_batch} lines the error is {error}");
            };
            assert_eq!((line, place), (5, VcfPlace::Line), "{lines_per_batch}");
            assert!(problem.contains("UTF-8"), "{lines_per_batch}: {problem}");
        }
    }

    #[test]
    fn an_alt_that_ends_in_a_comma_is_refused() {
        let vcf = vcf_of(&["chr1 100 . A T, . PASS . GT 0/0 0/1 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("ALT")));
    }

    #[test]
    fn a_ref_with_no_letter_in_it_is_refused() {
        // The REF of this line is empty: two tabs, one after the other.
        let vcf = vcf_of(&["chr1 100 .  T . PASS . GT 0/0 0/1 1/1"]);
        let error = error_reading(&vcf, VcfOptions::default());
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("REF")));
    }

    #[test]
    fn a_gzipped_source_cut_in_the_middle_of_a_member_is_refused() {
        let bytes = std::fs::read(reference_vcf("cases.vcf.gz")).unwrap();
        // The first member of this file is its header and the second its
        // four variants, so the cut is inside the second one. bgzip wrote
        // the file, so what the reader says of it is that it was cut short,
        // and not what the decoder found, an incomplete deflate stream.
        let cut = bytes.len().saturating_sub(20);
        let bytes = bytes.get(..cut).unwrap_or_default().to_vec();
        for lines_per_batch in BATCHES_TO_TRY {
            let mut reader =
                VcfReader::new(Cursor::new(bytes.clone()), VcfOptions::default()).unwrap();
            reader.set_lines_per_batch(lines_per_batch);
            assert_eq!(reader.individuals(), ["ind1", "ind2", "ind3"]);
            let error = rows_of(&mut reader).unwrap_err();
            assert!(
                matches!(error, Error::VcfBgzipEndMissing),
                "a file cut in the middle of a gzip member read in batches of \
                 {lines_per_batch} lines gives {error}"
            );
        }
    }

    #[test]
    fn a_gzip_that_bgzip_did_not_write_and_that_is_cut_is_an_error_of_the_input() {
        // Nothing says of such a file where it should end, so what a reader
        // has of it is what the decoder found: the bytes ran out.
        let plain = std::fs::read(reference_vcf("cases.vcf")).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, &plain).unwrap();
        let gzipped = encoder.finish().unwrap();
        let cut = gzipped.len().saturating_sub(20);
        let bytes = gzipped.get(..cut).unwrap_or_default().to_vec();

        let mut reader = VcfReader::new(Cursor::new(bytes), VcfOptions::default()).unwrap();
        let error = rows_of(&mut reader).unwrap_err();
        assert!(matches!(error, Error::Io(_)), "the error is {error}");
    }

    // A quality that is not finite, and a source that bgzip wrote and that
    // does not end with the mark of the end of a bgzipped file. Both are
    // errors that the owner decided on 20 September 2026, and pyNei reads
    // both files.

    #[test]
    fn a_quality_that_is_not_finite_is_refused() {
        // `1e400` is above what a float of 64 bits holds and `1e39` above
        // what one of 32 bits holds, which is what a block keeps, and both
        // read as an infinite quality. NaN is what a block holds for a
        // variant with no quality, so a NaN that was written in the file
        // would be read as a variant that has none.
        for quality in ["nan", "NaN", "inf", "-inf", "1e400", "1e39"] {
            let line = format!("chr1 100 . A T {quality} PASS . GT 0/0 0/1 1/1");
            let error = error_reading(&vcf_of(&[&line]), VcfOptions::default());
            let Error::VcfDataLine {
                line: number,
                place,
                problem,
            } = error
            else {
                panic!("the error of the quality `{quality}` is {error}");
            };
            assert_eq!(
                (number, place),
                (FIRST_DATA_LINE, VcfPlace::Column("QUAL")),
                "{quality}"
            );
            assert!(problem.contains(quality), "{quality}: {problem}");
        }
    }

    #[test]
    fn a_quality_that_is_not_finite_is_read_when_the_quality_was_not_asked_for() {
        // A column that is not parsed is not checked.
        let vcf = vcf_of(&["chr1 100 . A T nan PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::GTS | Needs::CHROM_POS);
        let rows = rows_of(&mut reader).expect("the rows");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].pos, 100);
        assert_eq!(rows[0].gts, [0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn a_quality_that_is_a_finite_number_is_read() {
        // The largest quality a block holds, and the smallest: a file of
        // real qualities has neither, and what they show is that the check
        // refuses the infinite ones and nothing else.
        let line = format!("chr1 100 . A T {} PASS . GT 0/0 0/1 1/1", f32::MAX);
        let rows = rows_read(&vcf_of(&[&line]), VcfOptions::default());
        assert_eq!(rows.first().map(|row| row.qual), Some(Some(f32::MAX)));
    }

    /// The bytes of one of the reference files, cut to `bytes` of them.
    fn cut_to(name: &str, bytes: usize) -> Vec<u8> {
        let whole = std::fs::read(reference_vcf(name)).unwrap();
        whole.get(..bytes).unwrap_or_default().to_vec()
    }

    /// The variants a reader gives before the error it ends with, which is
    /// what a source with no mark of its end gives: the variants first and
    /// the error where the reader would have said that there are no more.
    fn blocks_and_then_the_error(mut reader: VcfReader<Cursor<Vec<u8>>>) -> (Vec<usize>, Error) {
        let mut sizes = Vec::new();
        loop {
            match reader.next_block() {
                Ok(Some(block)) => sizes.push(block.num_vars),
                Ok(None) => panic!("the reader ended with no error after {sizes:?}"),
                Err(error) => return (sizes, error),
            }
        }
    }

    /// The last 28 bytes that went by to the decoder are what says whether
    /// a file that bgzip wrote ends with the mark of its end, and they are
    /// kept as the bytes are consumed, in chunks of whatever size the
    /// decoder asks for: fewer than 28, exactly 28, and more.
    #[test]
    fn the_last_bytes_of_a_source_are_the_last_of_what_went_by_in_any_chunks() {
        let bytes: Vec<u8> = (0..=255).collect();
        for chunk in [1, 27, 28, 40, 256] {
            let mut source = WithTheLastBytes::new(Cursor::new(bytes.clone()));
            let mut read = Vec::new();
            let mut left = bytes.len();
            while left > 0 {
                let take = chunk.min(left);
                let mut into = vec![0u8; take];
                std::io::Read::read_exact(&mut source, &mut into).expect("the bytes");
                read.extend_from_slice(&into);
                left = left.saturating_sub(take);
            }
            assert_eq!(read, bytes, "in chunks of {chunk}");
            // The last 28 of the 256 bytes, which is what a mark of the end
            // would be compared with.
            assert_eq!(source.last, bytes[228..], "in chunks of {chunk}");
        }

        // Fewer bytes than the mark have gone by: they are kept as they
        // are, and no mark is found in them.
        let mut source = WithTheLastBytes::new(Cursor::new(vec![7u8; 10]));
        let mut into = [0u8; 10];
        std::io::Read::read_exact(&mut source, &mut into).expect("the bytes");
        assert_eq!(source.last, [7u8; 10]);
        assert!(!source.ends_with_the_bgzip_mark());
    }

    #[test]
    fn a_bgzipped_source_without_the_mark_of_its_end_is_refused_after_its_variants() {
        // `many.vcf.gz` without the empty block of 28 bytes that bgzip
        // writes at the end of a file. Its 500 variants are given first,
        // in five blocks of 100, and the error comes where the reader
        // would have said that there are no more.
        let cut = cut_to("many.vcf.gz", 21904 - 28);
        let options = in_blocks_of(options(2, false), 100);
        let (sizes, error) =
            blocks_and_then_the_error(VcfReader::new(Cursor::new(cut), options).unwrap());
        assert_eq!(sizes, [100; 5]);
        let Error::VcfBgzipEndMissing = error else {
            panic!("the error is {error}");
        };
        let message = error.to_string();
        assert!(message.contains("bgzip"), "{message}");
        assert!(message.contains("cut"), "{message}");
    }

    #[test]
    fn a_bgzipped_source_cut_where_a_member_ends_gives_its_variants_and_then_the_error() {
        // The first 12336 bytes of `many.vcf.gz` are its two first gzip
        // members, 280 whole data lines: a decoder finds nothing wrong
        // there, and what says that the file is cut short is the mark of
        // the end that is not at its end.
        let cut = cut_to("many.vcf.gz", 12336);
        let options = in_blocks_of(options(2, false), 100);
        let (sizes, error) =
            blocks_and_then_the_error(VcfReader::new(Cursor::new(cut), options).unwrap());
        assert_eq!(sizes, [100, 100, 80]);
        assert!(
            matches!(error, Error::VcfBgzipEndMissing),
            "the error is {error}"
        );
    }

    /// A source that gives `bytes` bytes and then fails, which is the disc
    /// that a file is read from failing, and not a file that was cut short:
    /// the two are told apart by the error the source gives.
    struct FailsAfter<R: BufRead> {
        source: R,
        /// How many bytes are left before it fails.
        left: usize,
    }

    impl<R: BufRead> std::io::Read for FailsAfter<R> {
        fn read(&mut self, out: &mut [u8]) -> std::io::Result<usize> {
            let read = {
                let mut buffer = self.fill_buf()?;
                std::io::Read::read(&mut buffer, out)?
            };
            self.consume(read);
            Ok(read)
        }
    }

    impl<R: BufRead> BufRead for FailsAfter<R> {
        fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
            if self.left == 0 {
                return Err(std::io::Error::other("the disc of the test failed"));
            }
            let buffer = self.source.fill_buf()?;
            let take = buffer.len().min(self.left);
            Ok(&buffer[..take])
        }

        fn consume(&mut self, amount: usize) {
            self.left = self.left.saturating_sub(amount);
            self.source.consume(amount);
        }
    }

    #[test]
    fn a_bgzipped_source_cut_inside_a_member_is_the_error_of_the_mark_of_its_end() {
        // A cut inside a gzip member leaves the decoder without the bytes
        // it needs, where a cut where a member ends leaves it content: the
        // two are the same file for a user, one whose download stopped, and
        // the same error. The lines the decoder could give are given first.
        //
        // The members of `many.vcf.gz` end at the bytes 310, 12336, 21876
        // and 21904, so these three cuts are inside the second member, the
        // third, which is the last one with data lines, and the empty one
        // that marks the end.
        for (cut, blocks) in [
            (1000, vec![11]),
            (21000, vec![100, 100, 100, 100, 80]),
            (21890, vec![100; 5]),
        ] {
            let bytes = cut_to("many.vcf.gz", cut);
            let options = in_blocks_of(options(2, false), 100);
            let (sizes, error) =
                blocks_and_then_the_error(VcfReader::new(Cursor::new(bytes), options).unwrap());
            assert_eq!(sizes, blocks, "cut at {cut}");
            assert!(
                matches!(error, Error::VcfBgzipEndMissing),
                "cut at {cut}: the error is {error}"
            );
        }
    }

    #[test]
    fn a_source_that_fails_while_it_is_read_is_an_error_of_the_input() {
        // The disc that a file is read from failing is not a file that was
        // cut short: the blocks that were read are given and the error is
        // the one of the input, whatever the error of a source that ends
        // early would have been.
        let whole = std::fs::read(reference_vcf("many.vcf.gz")).unwrap();
        let source = FailsAfter {
            source: Cursor::new(whole),
            left: 13000,
        };
        let options = in_blocks_of(options(2, false), 100);
        let mut reader = VcfReader::new(source, options).unwrap();
        let mut sizes = Vec::new();
        let error = loop {
            match reader.next_block() {
                Ok(Some(block)) => sizes.push(block.num_vars),
                Ok(None) => panic!("the reader ended with no error after {sizes:?}"),
                Err(error) => break error,
            }
        };
        assert_eq!(sizes, [100, 100]);
        let Error::Io(error) = error else {
            panic!("the error is {error}");
        };
        assert!(error.to_string().contains("disc"), "{error}");
    }

    /// What says that bgzip wrote a source is the extra field `BC` in the
    /// header of its first gzip member, which is a flag that says that
    /// there is an extra field and the two bytes that name it. Neither half
    /// says it alone: a gzip with another extra field has the flag, and a
    /// gzip with no extra field can hold those two bytes where the name of
    /// one would be.
    #[test]
    fn a_gzip_with_another_extra_field_was_not_written_by_bgzip() {
        // A gzip header with the extra field `QQ` of two bytes: the ten
        // bytes of the header with the flag of an extra field, the length
        // of the field, its name, the length of its bytes and its bytes.
        let mut header = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff];
        header.extend_from_slice(&[0x06, 0x00, b'Q', b'Q', 0x02, 0x00, 0x00, 0x00]);
        assert!(!written_by_bgzip(&header), "the extra field `QQ`");

        // The same bytes with the flag taken off, which is a header that
        // holds no extra field at all whatever comes after it.
        let mut with_no_flag = header.clone();
        with_no_flag[GZIP_FLAGS] = 0x00;
        assert!(!written_by_bgzip(&with_no_flag), "no flag");

        // And the header of `many.vcf.gz`, which bgzip wrote.
        let bgzipped = std::fs::read(reference_vcf("many.vcf.gz")).unwrap();
        assert!(written_by_bgzip(&bgzipped), "many.vcf.gz");
        let mut without_the_flag = bgzipped.clone();
        without_the_flag[GZIP_FLAGS] = 0x00;
        assert!(!written_by_bgzip(&without_the_flag), "many.vcf.gz, no flag");
    }

    /// A gzip file whose first member carries an extra field that is not
    /// bgzip's is read to its end and asked for no mark: a reader that took
    /// the flag alone for bgzip's would refuse it.
    #[test]
    fn a_gzip_with_another_extra_field_is_read_to_its_end() {
        let plain = std::fs::read(reference_vcf("cases.vcf")).unwrap();
        let mut encoder = flate2::GzBuilder::new()
            .extra(vec![b'Q', b'Q', 0x02, 0x00, 0x00, 0x00])
            .write(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, &plain).unwrap();
        let gzipped = encoder.finish().unwrap();
        assert!(!written_by_bgzip(&gzipped));

        let mut reader = VcfReader::new(Cursor::new(gzipped), options(2, false)).unwrap();
        assert_eq!(rows_of(&mut reader).expect("the rows"), the_rows_of_cases());
    }

    #[test]
    fn a_gzipped_source_that_bgzip_did_not_write_is_read_to_its_end() {
        // A gzip file of one member, which has no mark of its end to
        // miss: the reader knows that a source was made by bgzip from the
        // extra field `BC` of its first member.
        let plain = std::fs::read(reference_vcf("many.vcf")).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, &plain).unwrap();
        let gzipped = encoder.finish().unwrap();
        assert!(gzipped.len() < plain.len());

        let mut reader = VcfReader::new(Cursor::new(gzipped), VcfOptions::default()).unwrap();
        let rows = rows_of(&mut reader).expect("the variants");
        assert_eq!(rows, rows_of_file("many.vcf", VcfOptions::default()));
        assert_eq!(rows.len(), 475);
    }

    // The blocks of a file: their sizes, and that what they hold does not
    // depend on the size, on the batches or on the threads.

    #[test]
    fn the_blocks_of_many_vcf_have_the_sizes_of_the_spec() {
        // "How it is verified" of `docs/specs/block.md`: five blocks of
        // 100, 100, 100, 100 and 75 with the default, five of 100 with
        // every variant given, and one of 475 with blocks of 1000.
        let sizes = |options: VcfOptions| -> Vec<usize> {
            let mut reader = reader_of_file("many.vcf", options);
            let blocks = blocks_of(&mut reader).expect("the blocks");
            blocks.iter().map(|block| block.num_vars).collect()
        };
        let default = VcfOptions::default();
        assert_eq!(sizes(in_blocks_of(default, 100)), [100, 100, 100, 100, 75]);
        assert_eq!(sizes(in_blocks_of(options(2, false), 100)), [100; 5]);
        assert_eq!(sizes(in_blocks_of(default, 1000)), [475]);
        // The 50 individuals of `many.vcf` and no size asked for: its 475
        // variants are one block, since the default is 10000 of them.
        assert_eq!(sizes(default), [475]);
    }

    #[test]
    fn the_blocks_of_many_vcf_hold_the_same_whatever_their_size_and_their_batches() {
        let expected = rows_of_file("many.vcf", VcfOptions::default());
        assert_eq!(expected.len(), 475);
        for num_vars_per_block in [1, 7, 100, 1000] {
            for lines_per_batch in [1, 3, 1024] {
                let mut reader = reader_of_file(
                    "many.vcf",
                    in_blocks_of(VcfOptions::default(), num_vars_per_block),
                );
                reader.set_lines_per_batch(lines_per_batch);
                let rows = rows_of(&mut reader).expect("the rows");
                assert_eq!(
                    rows, expected,
                    "blocks of {num_vars_per_block} variants read in batches of \
                     {lines_per_batch} lines"
                );
            }
        }
    }

    #[test]
    fn the_error_of_the_third_variant_comes_after_the_block_of_the_two_before_it() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T . PASS . GT 0/0 0/1 1/1",
            "chr1 300 rs3 A T . PASS . GT 0/0 0/1 0/0/1/1",
            "chr1 400 rs4 A T . PASS . GT 0/0 0/1 1",
        ]);
        let mut reader = reader_over(&vcf, in_blocks_of(VcfOptions::default(), 2));

        let first = reader.next_block().unwrap().expect("the first block");
        assert_eq!(first.num_vars, 2);
        assert_eq!(first.gts, [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]);

        // Two wrong lines, and the error is that of the first of them.
        let error = match reader.next_block() {
            Ok(block) => panic!("the reader gave {block:?}"),
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

        // The block that was being built is lost with the error, and there
        // is no block after it.
        assert!(reader.next_block().unwrap().is_none());
    }

    // The row parser: one data line, as bytes, into one row of a block.
    // What it gives has to be what the tables of "How it is verified" of
    // `docs/specs/io_vcf.md` have, and the errors of a data line the ones
    // that section lists. The cases that need a reader, the header, the
    // gzip, the lines that the FILTER leaves out and the blocks, are tested
    // where the reader is.

    /// The individuals of the header of the tests, by the names that the
    /// error of an individual carries.
    fn three_individuals() -> Vec<String> {
        ["ind1", "ind2", "ind3"]
            .iter()
            .map(|name| (*name).to_string())
            .collect()
    }

    /// A data line written with spaces where the file has tabs.
    fn data_line(line: &str) -> String {
        line.replace(' ', "\t")
    }

    /// What a slot of the genotypes of a row holds before the parse, which
    /// is no allele a VCF can give: a test that compares the genotypes sees
    /// a slot the parse did not write.
    const NOT_WRITTEN: i8 = i8::MIN;

    /// One data line parsed into a row, with its genotypes, as
    /// [`parse_row`] gives them.
    fn parse_the_bytes(
        line: &[u8],
        needs: Needs,
        ploidy: usize,
        individuals: &[String],
    ) -> Result<Row> {
        let rules = RowRules {
            needs,
            ploidy,
            individuals,
            panic_at_line: None,
        };
        let mut gts = if needs.contains(Needs::GTS) {
            vec![NOT_WRITTEN; individuals.len().saturating_mul(ploidy)]
        } else {
            Vec::new()
        };
        let mut row = ParsedRow::default();
        parse_row(line, FIRST_DATA_LINE, &rules, &mut gts, &mut row)?;
        Ok(Row {
            chrom: row.chrom.clone(),
            pos: row.pos,
            id: row.id.clone(),
            alleles: row.alleles().to_vec(),
            qual: if row.qual.is_nan() {
                None
            } else {
                Some(row.qual)
            },
            gts,
        })
    }

    /// One data line, written with spaces where the file has tabs, parsed
    /// into a row for the three individuals of the header of the tests.
    fn parse_the_line(line: &str, needs: Needs, ploidy: usize) -> Result<Row> {
        parse_the_bytes(
            data_line(line).as_bytes(),
            needs,
            ploidy,
            &three_individuals(),
        )
    }

    /// The row of a data line that the parser reads.
    fn row_of_the_line(line: &str, needs: Needs, ploidy: usize) -> Row {
        match parse_the_line(line, needs, ploidy) {
            Ok(row) => row,
            Err(error) => panic!("the line was refused: {error}"),
        }
    }

    /// The error of a data line that the parser refuses.
    fn error_of_the_line(line: &str, needs: Needs, ploidy: usize) -> Error {
        match parse_the_line(line, needs, ploidy) {
            Ok(row) => panic!("the line gave the row {row:?}"),
            Err(error) => error,
        }
    }

    #[test]
    fn the_four_lines_of_cases_vcf_are_parsed_into_their_rows() {
        let lines = [
            "chr1 100 rs1 A T 29.5 PASS . GT:DP 0/0:3 0/1:4 1/1:5",
            "chr1 200 . A T . q10 . GT:DP ./.:. 0|1:3 .|0:2",
            "chr1 300 . A G,T 67 PASS . GT 1/2 2|1 2/2",
            "chr1 400 . T . 47 PASS . GT 0/0 0/0 0/0",
        ];
        let rows: Vec<Row> = lines
            .iter()
            .map(|line| row_of_the_line(line, Needs::ALL, 2))
            .collect();
        assert_eq!(rows, the_rows_of_cases());
    }

    #[test]
    fn the_two_lines_of_differences_vcf_are_parsed_into_their_rows() {
        let lines = [
            "chr2 50 ms1 GTC G,GTCT 50 PASS . GT:DP 0/1:3 0/2 .",
            "chr2 60 . A <DEL>,* . PASS . GT /0/1 |2|2 0/0",
        ];
        let rows: Vec<Row> = lines
            .iter()
            .map(|line| row_of_the_line(line, Needs::ALL, 2))
            .collect();
        assert_eq!(rows, the_rows_of_differences());
    }

    #[test]
    fn a_genotype_of_another_ploidy_is_refused_with_its_individual() {
        for (line, genotype_ploidy) in [
            ("chr1 100 . A T . PASS . GT 0/0 0/0/1/1 1/1", 4),
            ("chr1 100 . A T . PASS . GT 0/0 1 1/1", 1),
        ] {
            let error = error_of_the_line(line, Needs::ALL, 2);
            let Error::VcfGenotypePloidy {
                line: number,
                individual,
                found,
                expected,
            } = error
            else {
                panic!("the error is {error}");
            };
            assert_eq!(
                (number, individual.as_str(), found, expected),
                (FIRST_DATA_LINE, "ind2", genotype_ploidy, 2)
            );
        }
    }

    #[test]
    fn a_tetraploid_line_read_with_the_ploidy_four_gives_four_alleles() {
        let row = row_of_the_line(
            "chr1 100 . A T . PASS . GT 0/0/1/1 0/1/1/1 ./././.",
            Needs::ALL,
            4,
        );
        assert_eq!(
            row.gts,
            vec![
                0,
                0,
                1,
                1,
                0,
                1,
                1,
                1,
                MISSING_ALLELE,
                MISSING_ALLELE,
                MISSING_ALLELE,
                MISSING_ALLELE
            ]
        );
    }

    #[test]
    fn an_allele_the_variant_does_not_declare_is_refused_with_the_genotypes_alone() {
        // The alleles of ALT are counted for every line, so that the allele
        // numbers are checked also when the texts of the alleles are not
        // kept: this line declares one alternative allele and its third
        // genotype carries the allele 2.
        let error = error_of_the_line("chr1 100 . A T . PASS . GT 0/0 0/1 1/2", Needs::GTS, 2);
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (line, place),
            (FIRST_DATA_LINE, VcfPlace::Individual("ind3".to_string()))
        );
    }

    #[test]
    fn the_largest_allele_popnei_holds_is_read() {
        // 127 alternative alleles and the reference one, and a genotype
        // that carries the last of them. 128 is the allele above it, which
        // the test below refuses.
        let mut alternatives = String::from("T");
        for _ in 1..127 {
            alternatives.push_str(",T");
        }
        let line = format!("chr1 100 . A {alternatives} . PASS . GT 0/0 127/0 1/1");
        let row = row_of_the_line(&line, Needs::ALL, 2);
        assert_eq!(row.alleles.len(), 128);
        assert_eq!(row.gts, [0, 0, 127, 0, 1, 1]);
    }

    #[test]
    fn an_allele_above_the_largest_one_is_refused() {
        let mut alternatives = String::from("T");
        for _ in 1..200 {
            alternatives.push_str(",T");
        }
        let line = format!("chr1 100 . A {alternatives} . PASS . GT 0/0 0/1 1/128");
        let error = error_of_the_line(&line, Needs::ALL, 2);
        let Error::VcfDataLine {
            line: number,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (number, place),
            (FIRST_DATA_LINE, VcfPlace::Individual("ind3".to_string()))
        );
        assert!(problem.contains("127"), "{problem}");
    }

    #[test]
    fn a_line_with_another_number_of_columns_of_individuals_is_refused() {
        for (line, count) in [
            ("chr1 100 . A T . PASS . GT 0/0 0/1", "2"),
            ("chr1 100 . A T . PASS . GT 0/0 0/1 1/1 0/0", "1"),
        ] {
            let error = error_of_the_line(line, Needs::ALL, 2);
            let Error::VcfDataLine {
                line: number,
                place,
                problem,
            } = error
            else {
                panic!("the error is {error}");
            };
            assert_eq!((number, place), (FIRST_DATA_LINE, VcfPlace::Line));
            assert!(problem.contains(count), "{problem}");
            assert!(problem.contains('3'), "{problem}");
        }
    }

    #[test]
    fn a_format_with_no_gt_is_refused_in_a_row() {
        let error = error_of_the_line("chr1 100 . A T . PASS . DP 3 4 5", Needs::ALL, 2);
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("FORMAT")));
    }

    #[test]
    fn a_column_of_an_individual_with_no_value_where_the_format_has_gt_is_refused() {
        let error = error_of_the_line("chr1 100 . A T . PASS . DP:GT 3:0/1 4 5:1/1", Needs::ALL, 2);
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (line, place),
            (FIRST_DATA_LINE, VcfPlace::Individual("ind2".to_string()))
        );
        // The column of that individual holds one value and the FORMAT
        // names two, so there is no genotype in it. A column whose value
        // where the GT is is empty is another line and another error, the
        // one of an allele number with no digit in it, which the test
        // below has.
        assert!(problem.contains("has no value"), "{problem}");
        assert!(problem.contains('4'), "{problem}");
    }

    #[test]
    fn a_column_of_an_individual_that_is_empty_is_not_an_allele_number() {
        // A line that ends in a tab has an empty column at its end, and a
        // line with two tabs one after another has one in the middle: both
        // are a column, and neither holds a genotype.
        for (line, individual) in [
            ("chr1 100 . A T . PASS . GT 0/0 0/1 ", "ind3"),
            ("chr1 100 . A T . PASS . GT 0/0  1/1", "ind2"),
        ] {
            let error = error_of_the_line(line, Needs::ALL, 2);
            let Error::VcfDataLine {
                line: number,
                place,
                problem,
            } = error
            else {
                panic!("the error of `{line}` is {error}");
            };
            assert_eq!(
                (number, place),
                (
                    FIRST_DATA_LINE,
                    VcfPlace::Individual(individual.to_string())
                ),
                "{line}"
            );
            assert!(problem.contains("is not an allele number"), "{problem}");
        }
    }

    #[test]
    fn a_position_that_is_not_a_number_is_refused_only_when_the_position_was_asked_for() {
        let line = "chr1 x . A T . PASS . GT 0/0 0/1 1/1";
        let error = error_of_the_line(line, Needs::ALL, 2);
        let Error::VcfDataLine {
            line: number,
            place,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((number, place), (FIRST_DATA_LINE, VcfPlace::Column("POS")));

        // A column that is not parsed is not checked.
        let row = row_of_the_line(line, Needs::GTS, 2);
        assert_eq!(row.pos, 0);
        assert_eq!(row.gts, vec![0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn a_quality_that_is_not_a_number_is_refused_only_when_the_quality_was_asked_for() {
        let line = "chr1 100 . A T x PASS . GT 0/0 0/1 1/1";
        let error = error_of_the_line(line, Needs::ALL, 2);
        let Error::VcfDataLine {
            line: number,
            place,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((number, place), (FIRST_DATA_LINE, VcfPlace::Column("QUAL")));

        let row = row_of_the_line(line, Needs::GTS, 2);
        assert_eq!(row.qual, None);
    }

    #[test]
    fn an_allele_with_no_letter_in_it_is_refused_in_ref_and_in_alt() {
        for (line, column) in [
            ("chr1 100 . A T, . PASS . GT 0/0 0/1 1/1", "ALT"),
            ("chr1 100 .  T . PASS . GT 0/0 0/1 1/1", "REF"),
        ] {
            let error = error_of_the_line(line, Needs::ALL, 2);
            let Error::VcfDataLine {
                line: number,
                place,
                ..
            } = error
            else {
                panic!("the error is {error}");
            };
            assert_eq!((number, place), (FIRST_DATA_LINE, VcfPlace::Column(column)));
        }
    }

    #[test]
    fn a_line_whose_bytes_are_not_valid_utf8_is_refused() {
        // The byte that is not text is in the REF column, as it is in the
        // test of the reader over the same case.
        let line = b"chr1\t200\t.\t\xffA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1";
        let error = match parse_the_bytes(line, Needs::ALL, 2, &three_individuals()) {
            Ok(row) => panic!("the line gave the row {row:?}"),
            Err(error) => error,
        };
        let Error::VcfDataLine {
            line: number,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((number, place), (FIRST_DATA_LINE, VcfPlace::Line));
        assert!(problem.contains("UTF-8"), "{problem}");
    }

    /// The columns of the individuals are read as bytes and are never text,
    /// which "Speed" of `docs/specs/io_vcf.md` asks for, so a byte that is
    /// not text in one of them is not found as bytes that are not UTF-8: it
    /// is a byte that is not a digit where an allele number is, and the
    /// error names the individual whose column it is in. The UTF-8 of the
    /// line is checked over the nine first columns, which are the ones
    /// whose text is kept.
    #[test]
    fn a_byte_that_is_not_text_in_the_column_of_an_individual_is_not_an_allele_number() {
        let line = b"chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/\xff\t1/1";
        let error = match parse_the_bytes(line, Needs::ALL, 2, &three_individuals()) {
            Ok(row) => panic!("the line gave the row {row:?}"),
            Err(error) => error,
        };
        let Error::VcfDataLine {
            line: number,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!(
            (number, place),
            (FIRST_DATA_LINE, VcfPlace::Individual("ind2".to_string()))
        );
        assert!(problem.contains("digits"), "{problem}");
        // The byte is shown as its number and not as the one character that
        // stands for every byte that could not be read: a user looks for it
        // in their file.
        assert!(problem.contains("\\xff"), "{problem}");
    }

    #[test]
    fn a_gt_that_is_not_the_first_key_of_the_format_is_read_into_a_row() {
        let row = row_of_the_line(
            "chr1 100 . A T . PASS . DP:GT 3:0/1 4:1/1 5:./.",
            Needs::ALL,
            2,
        );
        assert_eq!(row.gts, vec![0, 1, 1, 1, MISSING_ALLELE, MISSING_ALLELE]);
    }

    #[test]
    fn a_line_that_ends_in_an_end_of_line_is_read_as_one_that_does_not() {
        let line = data_line("chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1");
        let individuals = three_individuals();
        let bare = parse_the_bytes(line.as_bytes(), Needs::ALL, 2, &individuals).unwrap();
        for ending in ["\n", "\r\n"] {
            let ended = format!("{line}{ending}");
            let row = parse_the_bytes(ended.as_bytes(), Needs::ALL, 2, &individuals).unwrap();
            assert_eq!(row, bare, "the line that ends in {ending:?}");
        }
        assert_eq!(bare.gts, vec![0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn a_line_of_seven_columns_is_refused_when_the_genotypes_are_not_asked_for() {
        let error = error_of_the_line("chr1 100 . A T . PASS", Needs::ID | Needs::ALLELES, 2);
        let Error::VcfDataLine {
            line,
            place,
            problem,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
        assert!(problem.contains("INFO"), "{problem}");
    }

    #[test]
    fn a_format_of_dp_is_refused_when_the_genotypes_are_not_asked_for() {
        let error = error_of_the_line(
            "chr1 100 . A T . PASS . DP 3 4 5",
            Needs::ID | Needs::ALLELES,
            2,
        );
        let Error::VcfDataLine { line, place, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Column("FORMAT")));
    }

    #[test]
    fn a_line_with_no_column_of_an_individual_is_refused() {
        let error = error_of_the_line("chr1 100 . A T . PASS . GT", Needs::ID, 2);
        let Error::VcfDataLine {
            line,
            place,
            problem,
            ..
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (FIRST_DATA_LINE, VcfPlace::Line));
        assert!(problem.contains("individual"), "{problem}");
    }

    #[test]
    fn the_columns_of_the_individuals_are_not_read_when_the_genotypes_are_not_asked_for() {
        // A genotype of another ploidy and a line with the columns of two
        // individuals under a header with three: both are read, and `gts`
        // is empty.
        for line in [
            "chr1 100 rs1 A T . PASS . GT 0/0/1/1 0/1 1/1",
            "chr1 100 rs1 A T . PASS . GT 0/0 0/1",
        ] {
            let row = row_of_the_line(line, Needs::ID | Needs::ALLELES, 2);
            assert_eq!(row.id, "rs1");
            assert_eq!(row.alleles, ["A", "T"]);
            assert!(row.gts.is_empty(), "{line}");
        }
    }

    #[test]
    fn only_the_fields_that_were_asked_for_are_parsed() {
        let line = "chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1";

        let genotypes_alone = row_of_the_line(line, Needs::GTS, 2);
        assert_eq!(genotypes_alone.gts, vec![0, 0, 0, 1, 1, 1]);
        assert_eq!(genotypes_alone.chrom, "");
        assert_eq!(genotypes_alone.pos, 0);
        assert_eq!(genotypes_alone.id, "");
        assert!(genotypes_alone.alleles.is_empty());
        assert_eq!(genotypes_alone.qual, None);

        let texts = row_of_the_line(line, Needs::ID | Needs::ALLELES, 2);
        assert!(texts.gts.is_empty());
        assert_eq!(texts.id, "rs1");
        assert_eq!(texts.alleles, ["A", "T"]);
        assert_eq!(texts.chrom, "");
        assert_eq!(texts.qual, None);

        let places = row_of_the_line(line, Needs::CHROM_POS | Needs::QUAL, 2);
        assert_eq!(places.chrom, "chr1");
        assert_eq!(places.pos, 100);
        assert_eq!(places.qual, Some(29.5));
        assert_eq!(places.id, "");
    }

    #[test]
    fn an_allele_that_no_genotype_carries_is_read() {
        let row = row_of_the_line("chr1 100 . A G,T . PASS . GT 0/0 0/1 1/1", Needs::ALL, 2);
        assert_eq!(row.alleles, ["A", "G", "T"]);
        assert_eq!(row.gts, vec![0, 0, 0, 1, 1, 1]);
    }

    #[test]
    fn a_genotype_written_as_a_dot_is_missing_in_every_allele_of_the_ploidy() {
        let row = row_of_the_line(
            "chr1 100 . A T . PASS . GT . 0/1/1/1 ./././.",
            Needs::GTS,
            4,
        );
        assert_eq!(row.gts.get(..4), Some([MISSING_ALLELE; 4].as_slice()));
        // A genotype of a dot for each allele is the same genotype written
        // out, and one of two dots under the ploidy 4 is of another ploidy.
        let error = error_of_the_line("chr1 100 . A T . PASS . GT . ./. 0/0/0/0", Needs::GTS, 4);
        let Error::VcfGenotypePloidy { individual, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(individual, "ind2");
    }

    /// The row the parser is given holds one allele for each individual of
    /// the header times the ploidy. A caller that gives it another number
    /// has a defect, and its genotypes would be read one at the place of
    /// another, so the parse refuses it instead of filling what it can.
    #[test]
    fn a_row_of_genotypes_that_is_not_of_the_size_of_the_line_is_refused() {
        let individuals = three_individuals();
        let rules = RowRules {
            needs: Needs::GTS,
            ploidy: 2,
            individuals: &individuals,
            panic_at_line: None,
        };
        let line = data_line("chr1 100 . A T . PASS . GT 0/0 0/1 1/1");
        let mut gts = vec![NOT_WRITTEN; 4];
        let mut row = ParsedRow::default();
        let error =
            parse_row(line.as_bytes(), FIRST_DATA_LINE, &rules, &mut gts, &mut row).unwrap_err();
        let Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((array, found, expected), ("gts", 4, 6));
    }

    // `many.vcf`, against what bcftools 1.24 read in it and the counts of
    // "How it is verified" of `docs/specs/io_vcf.md`.

    /// The ploidy of the 50 individuals of `many.vcf`, which
    /// `tests/reference/vcf/make_reference.py` writes as diploid.
    const MANY_PLOIDY: usize = 2;

    /// How many columns the output of `bcftools query` has before the
    /// genotypes: CHROM, POS, ID, REF, ALT, QUAL and FILTER, the seven of
    /// the format that `make_reference.py` gave bcftools.
    const COLUMNS_BEFORE_THE_GENOTYPES: usize = 7;

    /// Where the FILTER is among them, counted from 0.
    const FILTER_COLUMN: usize = 6;

    /// What the tests of `many.vcf` compare, one variant of it: the name of
    /// its chromosome, its position and the alleles of every genotype, the
    /// ploidy of them for each individual.
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct Site {
        chrom: String,
        pos: u64,
        gts: Vec<i8>,
    }

    /// One line of `many.bcftools.tsv`, what bcftools 1.24 printed for one
    /// variant of `many.vcf`.
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct ReferenceRow {
        site: Site,
        /// Whether its FILTER column is `PASS` or a dot, which is what the
        /// reader gives by default.
        passed: bool,
    }

    /// The alleles of a genotype as bcftools prints it, `0|1` or `./.`: the
    /// numbers of the alleles the variant declares, and [`MISSING_ALLELE`]
    /// for the dot of an allele that was not called.
    fn alleles_of(genotype: &str) -> impl Iterator<Item = i8> + '_ {
        genotype.split(['/', '|']).map(|allele| {
            if allele == MISSING_VALUE {
                MISSING_ALLELE
            } else {
                allele.parse().unwrap()
            }
        })
    }

    /// The rows of the file that `bcftools query` wrote beside a reference
    /// VCF, one line per variant with the seven columns above and then the
    /// GT of every individual.
    fn reference_rows(name: &str) -> Vec<ReferenceRow> {
        let text = std::fs::read_to_string(reference_vcf(name))
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let mut rows = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            let columns: Vec<&str> = line.split('\t').collect();
            assert!(
                columns.len() > COLUMNS_BEFORE_THE_GENOTYPES,
                "{name}: `{line}` has {} columns",
                columns.len()
            );
            rows.push(ReferenceRow {
                site: Site {
                    chrom: columns[0].to_string(),
                    pos: columns[1].parse().unwrap(),
                    gts: columns[COLUMNS_BEFORE_THE_GENOTYPES..]
                        .iter()
                        .flat_map(|genotype| alleles_of(genotype))
                        .collect(),
                },
                passed: matches!(columns[FILTER_COLUMN], "PASS" | MISSING_VALUE),
            });
        }
        rows
    }

    /// The chromosome, the position and the genotypes of the variants the
    /// reader gave, with the name of the chromosome that `rows_of` already
    /// took from the table of the reader.
    fn sites_of(rows: &[Row]) -> Vec<Site> {
        rows.iter()
            .map(|row| Site {
                chrom: row.chrom.clone(),
                pos: row.pos,
                gts: row.gts.clone(),
            })
            .collect()
    }

    /// The rows of bcftools the reader has to give with `options`: all of
    /// them, or the ones whose FILTER passed.
    fn expected_sites(rows: &[ReferenceRow], options: VcfOptions) -> Vec<Site> {
        rows.iter()
            .filter(|row| !options.only_passed || row.passed)
            .map(|row| row.site.clone())
            .collect()
    }

    /// The variants one by one, so that a file of 500 variants that differs
    /// in one of them says which one and not that two long lists differ.
    fn assert_the_same_sites(given: &[Site], expected: &[Site], read: &str) {
        assert_eq!(
            given.len(),
            expected.len(),
            "{read}: the number of variants"
        );
        for (index, (given, expected)) in given.iter().zip(expected).enumerate() {
            assert_eq!(
                given, expected,
                "{read}: the variant {index}, counted from 0"
            );
        }
    }

    /// The eight counts of the table of `many.vcf` of "How it is verified"
    /// of `docs/specs/io_vcf.md`, counted in the variants the reader gave.
    #[derive(Debug, PartialEq, Eq)]
    struct Counts {
        variants: u64,
        variants_in_chr2: u64,
        variants_with_two_alternative_alleles: u64,
        /// The genotypes with an allele that was not called, the half
        /// called ones among them.
        missing_genotypes: u64,
        /// The genotypes with an allele called and another one not.
        half_called_genotypes: u64,
        missing_alleles: u64,
        called_alleles: u64,
        /// The sum of the numbers of the alleles that were called, 0 for
        /// the reference allele and 1 and 2 for the alternative ones. It is
        /// the count that a genotype read at the wrong allele changes.
        sum_of_the_called_alleles: u64,
    }

    /// A count of the tests as the number the table of the spec has. Every
    /// count of a file of 500 variants of 50 individuals fits in a `u64`.
    fn counted(number: usize) -> u64 {
        u64::try_from(number).unwrap()
    }

    fn counts_of(rows: &[Row], ploidy: usize) -> Counts {
        let genotypes = || rows.iter().flat_map(|row| row.gts.chunks_exact(ploidy));
        let alleles = || rows.iter().flat_map(|row| row.gts.iter());
        let called = |allele: &&i8| **allele != MISSING_ALLELE;
        let missing_in = |genotype: &[i8]| {
            genotype
                .iter()
                .filter(|allele| **allele == MISSING_ALLELE)
                .count()
        };
        Counts {
            variants: counted(rows.len()),
            variants_in_chr2: counted(rows.iter().filter(|row| row.chrom == "chr2").count()),
            // The reference allele and the two alternative ones.
            variants_with_two_alternative_alleles: counted(
                rows.iter().filter(|row| row.alleles.len() == 3).count(),
            ),
            missing_genotypes: counted(
                genotypes()
                    .filter(|genotype| missing_in(genotype) != 0)
                    .count(),
            ),
            half_called_genotypes: counted(
                genotypes()
                    .filter(|genotype| {
                        let missing = missing_in(genotype);
                        missing != 0 && missing != genotype.len()
                    })
                    .count(),
            ),
            missing_alleles: counted(
                alleles()
                    .filter(|allele| **allele == MISSING_ALLELE)
                    .count(),
            ),
            called_alleles: counted(alleles().filter(called).count()),
            sum_of_the_called_alleles: alleles()
                .filter(called)
                .map(|allele| u64::from(allele.unsigned_abs()))
                .sum(),
        }
    }

    #[test]
    fn every_variant_of_many_vcf_is_read_as_bcftools_read_it() {
        let reference = reference_rows("many.bcftools.tsv");
        assert_eq!(reference.len(), 500);
        let options = options(MANY_PLOIDY, false);
        let expected = expected_sites(&reference, options);
        for name in ["many.vcf", "many.vcf.gz"] {
            let given = sites_of(&rows_of_file(name, options));
            assert_the_same_sites(&given, &expected, &format!("{name}, every variant"));
        }
    }

    #[test]
    fn the_default_gives_the_variants_of_many_vcf_whose_filter_passed() {
        let reference = reference_rows("many.bcftools.tsv");
        let expected = expected_sites(&reference, VcfOptions::default());
        // `bcftools view -f .,PASS` left 475 of the 500 variants.
        assert_eq!(expected.len(), 475);
        for name in ["many.vcf", "many.vcf.gz"] {
            let given = sites_of(&rows_of_file(name, VcfOptions::default()));
            assert_the_same_sites(&given, &expected, &format!("{name}, the default"));
        }
    }

    #[test]
    fn the_first_genotypes_of_many_vcf_are_the_ones_of_the_spec() {
        // The spec gives the first five genotypes of the first variant,
        // chr1 1000, as `1/1 .|. 1/0 1/1 0/1`.
        let first_five = [1, 1, MISSING_ALLELE, MISSING_ALLELE, 1, 0, 1, 1, 0, 1];
        for name in ["many.vcf", "many.vcf.gz"] {
            let rows = rows_of_file(name, VcfOptions::default());
            let first = &rows[0];
            assert_eq!((first.chrom.as_str(), first.pos), ("chr1", 1000), "{name}");
            assert_eq!(
                first.gts.get(..first_five.len()),
                Some(&first_five[..]),
                "{name}"
            );
        }
    }

    #[test]
    fn the_counts_of_many_vcf_with_every_variant_given_are_the_ones_of_the_spec() {
        let expected = Counts {
            variants: 500,
            variants_in_chr2: 250,
            variants_with_two_alternative_alleles: 54,
            missing_genotypes: 1511,
            half_called_genotypes: 257,
            missing_alleles: 2765,
            called_alleles: 47235,
            sum_of_the_called_alleles: 25954,
        };
        for name in ["many.vcf", "many.vcf.gz"] {
            let rows = rows_of_file(name, options(MANY_PLOIDY, false));
            assert_eq!(counts_of(&rows, MANY_PLOIDY), expected, "{name}");
        }
    }

    #[test]
    fn the_counts_of_many_vcf_with_the_default_are_the_ones_of_the_spec() {
        let expected = Counts {
            variants: 475,
            variants_in_chr2: 238,
            variants_with_two_alternative_alleles: 53,
            missing_genotypes: 1431,
            half_called_genotypes: 240,
            missing_alleles: 2622,
            called_alleles: 44878,
            sum_of_the_called_alleles: 24831,
        };
        for name in ["many.vcf", "many.vcf.gz"] {
            let rows = rows_of_file(name, VcfOptions::default());
            assert_eq!(counts_of(&rows, MANY_PLOIDY), expected, "{name}");
        }
    }

    /// The variants of `many.vcf`, each with the number its reader gave to
    /// its chromosome, read inside a rayon pool of `threads` threads, in
    /// blocks of `num_vars_per_block` variants and in batches of
    /// `lines_per_batch` lines.
    ///
    /// The pool is built here and is not rayon's global one, which has one
    /// thread per core of the machine: `install` runs the reader on this
    /// one instead. rayon is a dependency of the targets that are not wasm,
    /// so this and what uses it are compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    fn many_vcf_read_in_a_pool(
        threads: usize,
        num_vars_per_block: usize,
        lines_per_batch: usize,
    ) -> Vec<(u32, Row)> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut reader = reader_of_file(
                "many.vcf",
                in_blocks_of(options(MANY_PLOIDY, false), num_vars_per_block),
            );
            reader.set_lines_per_batch(lines_per_batch);
            let mut variants = Vec::new();
            loop {
                let block = match reader.next_block() {
                    Ok(Some(block)) => block,
                    Ok(None) => return variants,
                    Err(error) => panic!("many.vcf: the reader stopped at {error}"),
                };
                for view in block.variants() {
                    let number = view.chrom().expect("the chromosome");
                    variants.push((
                        number,
                        Row {
                            chrom: reader.chroms().name(number).unwrap_or_default().to_string(),
                            pos: view.pos().expect("the position"),
                            id: view.id().expect("the id").to_string(),
                            alleles: (0..view.num_alleles().expect("the alleles"))
                                .map(|allele| view.allele(allele).unwrap_or_default().to_string())
                                .collect(),
                            qual: view.qual().filter(|qual| !qual.is_nan()),
                            gts: view.gts().to_vec(),
                        },
                    ));
                }
            }
        })
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn many_vcf_gives_the_same_variants_in_a_pool_of_one_thread_and_in_one_of_four() {
        // 64 lines a batch and blocks of 100 variants, and the file has 500
        // data lines, so the pool of four threads parses several batches
        // and every thread parses lines of each.
        let one_thread = many_vcf_read_in_a_pool(1, 100, 64);
        let four_threads = many_vcf_read_in_a_pool(4, 100, 64);
        assert_eq!(one_thread.len(), 500);
        assert_eq!(one_thread, four_threads);

        // They are the variants bcftools read, and the numbers of the
        // chromosomes are the ones of the order of the file, whose first
        // name is `chr1`.
        let rows: Vec<Row> = four_threads.iter().map(|(_, row)| row.clone()).collect();
        let reference = reference_rows("many.bcftools.tsv");
        let expected = expected_sites(&reference, options(MANY_PLOIDY, false));
        assert_the_same_sites(&sites_of(&rows), &expected, "many.vcf in a pool of four");
        for (number, row) in &four_threads {
            let expected_number = if row.chrom == "chr1" { 0 } else { 1 };
            assert_eq!(*number, expected_number, "{} {}", row.chrom, row.pos);
        }
    }

    #[test]
    fn a_reader_whose_parse_did_not_come_back_gives_an_error_and_no_block() {
        // Eight data lines, in blocks of four, and the parse of the second
        // line of the second block panics: the first block is given, and
        // then the call that reads the second panics inside the parse and
        // unwinds through the reader.
        let lines: Vec<String> = (1..=8)
            .map(|pos| format!("chr1 {pos}00 . A T . PASS . GT 0/0 0/1 1/1"))
            .collect();
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let mut reader = reader_over(&vcf_of(&lines), in_blocks_of(VcfOptions::default(), 4));
        // The three lines of the header come first, so the sixth data line
        // is the line 9 of the file.
        reader.panic_at_line(9);

        let mut blocks = 0;
        let hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let went_on = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            while reader.next_block().unwrap().is_some() {
                blocks += 1;
            }
        }));
        std::panic::set_hook(hook);
        assert!(went_on.is_err(), "the parse did not panic");
        assert_eq!(blocks, 1);

        // The lines of the batch that was being parsed were never parsed,
        // and a reader that went on would give the variants of the ones
        // that were and drop the others without a word.
        let error = reader.next_block().unwrap_err();
        let Error::VcfParseNotFinished { line } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(line, 11);
        // Every call after it gives no block and not an error again: a
        // reader ends at its error, and the one that was given is the one
        // that says what happened.
        assert!(reader.next_block().unwrap().is_none());
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn a_wrong_line_of_a_later_batch_comes_after_the_blocks_that_were_read_before_it() {
        let mut lines: Vec<String> = (1..=200)
            .map(|pos| format!("chr1 {pos} . A T . PASS . GT 0/0 0/1 1/1"))
            .collect();
        // A haploid genotype in a file read as diploid, in the 201st data
        // line, which is the line 204 of the file and, with batches of 64
        // lines, is in the fourth batch and not in the first.
        lines.push("chr1 201 . A T . PASS . GT 1 0 1".to_string());
        lines.push("chr1 202 . A T . PASS . GT 0/0 0/1 1/1".to_string());
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let vcf = vcf_of(&lines);

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            let mut reader = reader_over(&vcf, in_blocks_of(VcfOptions::default(), 100));
            reader.set_lines_per_batch(64);
            for block in 1..=2 {
                let given = reader
                    .next_block()
                    .unwrap_or_else(|error| panic!("the block {block}: {error}"))
                    .expect("a block");
                assert_eq!(given.num_vars, 100, "the block {block}");
            }

            let error = reader.next_block().unwrap_err();
            let Error::VcfGenotypePloidy {
                line,
                individual,
                found,
                expected,
            } = error
            else {
                panic!("the error is {error}");
            };
            assert_eq!(
                (line, individual.as_str(), found, expected),
                (204, "ind1", 1, 2)
            );
            // The good line after the wrong one was parsed in the same
            // batch and is not given: an error ends the reader.
            assert!(reader.next_block().unwrap().is_none());
        });
    }

    /// The data lines of a reference VCF as a batch, with the bytes of the
    /// lines and the rows they are parsed into, which is what the reader
    /// gives the parse.
    fn a_batch_of(name: &str) -> (Vec<u8>, Vec<BatchRow>) {
        let text = std::fs::read_to_string(reference_vcf(name))
            .unwrap_or_else(|error| panic!("{name}: {error}"));
        let mut bytes = Vec::new();
        let mut batch = Vec::new();
        for (index, line) in text.lines().enumerate() {
            if line.starts_with('#') {
                continue;
            }
            let start = bytes.len();
            bytes.extend_from_slice(line.as_bytes());
            let mut row = BatchRow::new();
            row.line = start..bytes.len();
            // The lines of a file are counted from 1.
            row.number = u64::try_from(index).unwrap().saturating_add(1);
            batch.push(row);
        }
        (bytes, batch)
    }

    /// What one line of a batch was parsed into, to be compared with what
    /// the other way of parsing gave for the same line.
    fn parsed(
        row: &BatchRow,
        gts: &[i8],
    ) -> (String, String, u64, Vec<i8>, String, Vec<String>, u32) {
        (
            row.error
                .as_ref()
                .map_or_else(|| "no error".to_string(), Error::to_string),
            row.row.chrom.clone(),
            row.row.pos,
            gts.to_vec(),
            row.row.id.clone(),
            row.row.alleles().to_vec(),
            row.row.qual.to_bits(),
        )
    }

    #[test]
    fn the_lines_parsed_one_after_another_give_what_the_threads_give() {
        let individuals = reader_of_file("many.vcf", VcfOptions::default())
            .individuals()
            .to_vec();
        let rules = RowRules {
            needs: Needs::ALL,
            ploidy: MANY_PLOIDY,
            individuals: &individuals,
            panic_at_line: None,
        };
        let gts_per_variant = individuals.len().saturating_mul(MANY_PLOIDY);

        let (text, mut on_the_threads) = a_batch_of("many.vcf");
        let (_, mut one_by_one) = a_batch_of("many.vcf");
        assert_eq!(on_the_threads.len(), 500);
        let mut gts_of_the_threads = vec![MISSING_ALLELE; on_the_threads.len() * gts_per_variant];
        let mut gts_one_by_one = gts_of_the_threads.clone();

        parse_rows(&mut on_the_threads, &text, &mut gts_of_the_threads, &rules);
        parse_rows_one_by_one(&mut one_by_one, &text, &mut gts_one_by_one, &rules);

        for (index, (threads, one)) in on_the_threads.iter().zip(&one_by_one).enumerate() {
            let row = index * gts_per_variant..(index + 1) * gts_per_variant;
            assert_eq!(
                parsed(threads, &gts_of_the_threads[row.clone()]),
                parsed(one, &gts_one_by_one[row]),
                "the line {index}, counted from 0"
            );
        }
        let wrong = one_by_one.iter().filter(|row| row.error.is_some()).count();
        assert_eq!(wrong, 0);
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn a_batch_holds_more_than_one_line_where_there_are_threads() {
        // Every test of this file passes with one line in a batch: what a
        // block holds does not depend on how many lines were parsed
        // together or on how many threads parsed them, which is what the
        // reader promises. So nothing here notices a reader that stopped
        // reading ahead, and what would notice is the time it takes: the
        // benchmark `benches/read_vcf.rs` on the 400 MB VCF of 100000
        // variants of 1000 individuals takes about a second on one thread
        // and a sixth of that on the 18 cores of the owner's machine, and
        // with one line in a batch there is nothing for the other 17 to do.
        // The assertion is over a constant, and clippy asks for a const
        // block, which would refuse to compile instead of failing as a
        // test: a test that fails is what names this file and this reason.
        let lines_per_batch = std::hint::black_box(LINES_PER_BATCH);
        assert!(
            lines_per_batch > 1,
            "a batch of {lines_per_batch} line gives the threads of rayon one line to share"
        );
    }

    #[test]
    fn a_bound_of_bytes_smaller_than_the_file_cuts_the_batches_and_changes_no_result() {
        let expected = rows_of_file("many.vcf", VcfOptions::default());
        assert_eq!(expected.len(), 475);
        // `many.vcf` has 500 data lines of about 230 bytes each, read here
        // in one block of 1000 variants. A bound of one byte holds one line
        // in every batch, which is the line a batch holds whatever the
        // bound says; the batch of the last line stops at the bound without
        // seeing the end of the file, so one more batch is filled, which
        // reads no line: 501. The bound of the code, 8 MiB, takes the 500
        // lines and the end of the file in one.
        for (bytes, batches) in [(1, 501), (BYTES_PER_BATCH, 1)] {
            let mut reader = reader_of_file("many.vcf", in_blocks_of(VcfOptions::default(), 1000));
            reader.set_bytes_per_batch(bytes);
            let rows = rows_of(&mut reader).unwrap();
            assert_eq!(rows, expected, "with {bytes} bytes a batch");
            assert_eq!(
                reader.batches_filled(),
                batches,
                "with {bytes} bytes a batch"
            );
        }
    }

    #[test]
    fn a_file_read_one_line_at_a_time_gives_what_it_gives_in_one_batch() {
        // A batch of one line is what wasm reads, where there are no
        // threads, and the cargo tests are run natively.
        let mut reader = reader_of_file("many.vcf", VcfOptions::default());
        reader.set_lines_per_batch(1);
        let one_line_at_a_time = rows_of(&mut reader).unwrap();
        assert_eq!(
            one_line_at_a_time,
            rows_of_file("many.vcf", VcfOptions::default())
        );
        assert_eq!(one_line_at_a_time.len(), 475);
    }

    /// The rows of a batch are the buffers that the lines of the next batch
    /// are parsed into, so a reader of a file of any length keeps the rows
    /// of one batch and no more: the reader that this one took the place of
    /// kept the buffers of one variant for each line of a batch, and the
    /// test of that was over the buffer of the genotypes of a variant.
    #[test]
    fn the_rows_of_a_batch_are_the_ones_the_next_batch_is_parsed_into() {
        let mut reader = reader_of_file("many.vcf", in_blocks_of(VcfOptions::default(), 1000));
        reader.set_lines_per_batch(8);
        let rows = rows_of(&mut reader).unwrap();
        assert_eq!(rows.len(), 475);
        assert_eq!(reader.batch.len(), 8);
        assert!(reader.batches_filled() > 60, "{}", reader.batches_filled());
        // The text of the lines of a batch is written over by the next
        // batch, so what the reader holds is the lines of one batch and not
        // the file: `many.vcf` is 110879 bytes and eight of its lines are
        // about 1900.
        assert!(
            reader.text.capacity() < 8192,
            "the text of a batch of eight lines holds {} bytes",
            reader.text.capacity()
        );
        // And what bounds it for a file of many individuals, whose lines
        // are longer, is the bound in bytes.
        let mut reader = reader_of_file("many.vcf", in_blocks_of(VcfOptions::default(), 1000));
        reader.set_bytes_per_batch(1024);
        assert_eq!(rows_of(&mut reader).unwrap().len(), 475);
        assert!(
            reader.text.capacity() < 4096,
            "with a bound of 1024 bytes a batch holds {} bytes",
            reader.text.capacity()
        );
    }
}
