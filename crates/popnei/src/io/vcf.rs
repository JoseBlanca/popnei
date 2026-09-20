//! The VCF reader: it gives the variants of a VCF, plain or gzipped, one at
//! a time.
//!
//! [`VcfReader::new`] takes any source of bytes and reads the header of the
//! VCF, so the individuals are known before a variant is, and
//! [`VcfReader::from_path`] does the same for a caller that has a path. The
//! reader is a [`VariantReader`]: the consumer owns one [`Variant`] and
//! lends it to `read_variant` again and again.
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
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

use crate::error::{Error, Result};
use crate::variant::{ChromTable, MAX_ALLELE, MISSING_ALLELE, Needs, Variant, VariantReader};

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
/// Nobody has measured it. What bounds it from above is the memory a batch
/// holds, the text of its lines and the variants they were parsed into: on
/// the VCF of "Speed" of `docs/specs/io_vcf.md`, 100000 variants of 1000
/// individuals in 400 MB, a line is 4 KB and its genotypes 2000 alleles, so
/// 1024 lines are 4 MB of text and 2 MB of genotypes. What bounds it from
/// below is that the lines of one batch are what the threads share: on the
/// 18 cores of the machine of that section, 1024 lines are 56 lines a
/// thread, so a thread that starts late costs the others a fraction of the
/// batch and not a wait for half of it. The measurement of the reader on
/// that file is task 5.2 of `docs/plans/vcf-to-blocks.md`, and it is what
/// says whether this is the right number.
#[cfg(not(target_family = "wasm"))]
const LINES_PER_BATCH: usize = 1024;

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
/// in lines alone lets the memory of a reader grow with the panel: with
/// the 1024 lines of [`LINES_PER_BATCH`], a reader of the VCF of "Speed"
/// of `docs/specs/io_vcf.md`, 1000 individuals, holds 13.3 MB of resident
/// memory, one of 10000 individuals 82.7 MB and one of 100000 individuals
/// near 0.8 GB, which the review of work package 5 of
/// `docs/plans/vcf-to-blocks.md` measured. The text of the lines is what
/// that memory is made of, the variants parsed from them and the growth of
/// the buffers by doubling, so bounding the text bounds all of it, at
/// about three times this number.
///
/// 8 MiB is twice the 4.1 MB that 1024 lines of that 400 MB file hold, so
/// the batches of a file of a thousand individuals are the 1024 lines they
/// were and the timings of task 5.2 hold; a file of 10000 individuals gets
/// about 200 lines in a batch and one of 100000 about 20. A batch holds
/// one line whatever its bytes are.
const BYTES_PER_BATCH: usize = 8 * 1024 * 1024;

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
/// of the file has, and whether the variants that failed a filter are left
/// out.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfOptions {
    /// How many alleles every genotype has. 1 or more.
    pub ploidy: usize,
    /// Skip the variants whose FILTER is neither PASS nor a dot.
    pub only_passed: bool,
}

impl Default for VcfOptions {
    /// [`DEFAULT_PLOIDY`] and [`DEFAULT_ONLY_PASSED`], which are 2 and the
    /// variants that passed their filters alone.
    fn default() -> VcfOptions {
        VcfOptions {
            ploidy: DEFAULT_PLOIDY,
            only_passed: DEFAULT_ONLY_PASSED,
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

/// The bytes of the VCF, decompressed when the source was gzipped.
enum VcfSource<R: BufRead> {
    /// The source as it was given, with its first bytes in front of it.
    Plain(WithFirstBytes<R>),
    /// The source through the decoder of gzip, which goes on to the next
    /// member of the file when one ends.
    Gzipped(BufReader<MultiGzDecoder<WithFirstBytes<R>>>),
}

impl<R: BufRead> VcfSource<R> {
    /// The bytes that are already in the buffer, without consuming them.
    fn first_bytes(&mut self) -> std::io::Result<&[u8]> {
        match self {
            VcfSource::Plain(source) => source.fill_buf(),
            VcfSource::Gzipped(source) => source.fill_buf(),
        }
    }

    /// The next line, with its end of line, appended to `line`. 0 at the
    /// end of the source.
    ///
    /// Bytes that are not valid UTF-8 come back as an error of the input
    /// whose kind says so, which the callers turn into an error of the line
    /// they were reading, since a VCF is text and the number of the line is
    /// what a user needs.
    fn read_line(&mut self, line: &mut String) -> std::io::Result<usize> {
        match self {
            VcfSource::Plain(source) => source.read_line(line),
            VcfSource::Gzipped(source) => source.read_line(line),
        }
    }
}

/// What one data line of the batch gave when it was parsed.
enum LineOutcome {
    /// No variant: the line was empty, or its FILTER says that the variant
    /// failed a filter and the options leave those out. It is also what a
    /// line that was handed out is left with.
    NoVariant,
    /// The variant in the `var` of the line.
    Variant,
    /// The line is wrong, and this is the error that the reader gives at
    /// the read that would have given its variant.
    Wrong(Error),
}

/// One line of the batch the reader parses, with the variant it gave and
/// the buffers that the next line read into this place writes over.
///
/// Every line of a batch has its own, so that the threads that parse them
/// share nothing.
struct BatchLine {
    /// The line as the source gave it, with its end of line.
    text: String,
    /// Its number in the file, counted from 1 with the lines of the header.
    number: u64,
    /// The name of the chromosome of its variant, which gets its number
    /// when the variant is handed out.
    chrom_name: String,
    /// The variant the line was parsed into, whose `chrom` is not given
    /// yet.
    var: Variant,
    /// The strings of the alleles of the variants parsed here before,
    /// written over by the next one that asks for the alleles.
    spare_alleles: Vec<String>,
    /// What the line gave.
    outcome: LineOutcome,
}

impl BatchLine {
    /// A line with no text and no variant, which the reader adds to the
    /// batch when it reads more lines than the batch has held so far.
    fn new() -> BatchLine {
        BatchLine {
            text: String::new(),
            number: 0,
            chrom_name: String::new(),
            var: Variant::new(),
            spare_alleles: Vec::new(),
            outcome: LineOutcome::NoVariant,
        }
    }

    /// The variant of the text of this line, into its own variant, and what
    /// it gave into its `outcome`.
    ///
    /// The variant is cleared first and the strings of the alleles it held
    /// are kept, so that a line parsed again, and every line of the batches
    /// that come after this one, allocates nothing.
    fn parse(&mut self, rules: &ParseRules<'_>) {
        let BatchLine {
            text,
            number,
            chrom_name,
            var,
            spare_alleles,
            outcome,
        } = self;
        #[cfg(test)]
        tests::panic_if_the_test_asked_for_it(*number, rules);
        var.clear_but_the_alleles();
        spare_alleles.append(&mut var.alleles);
        chrom_name.clear();
        *outcome = match parse_data_line(
            without_the_line_end(text),
            *number,
            rules,
            var,
            chrom_name,
            spare_alleles,
        ) {
            Ok(true) => LineOutcome::Variant,
            Ok(false) => LineOutcome::NoVariant,
            Err(error) => LineOutcome::Wrong(error),
        };
    }
}

/// The lines of a batch, each parsed into the variant of its own line.
///
/// Natively they are parsed on the threads of rayon's global pool: no two
/// lines share a buffer, and what depends on the order of the file, the
/// number of the chromosome, is given later, when the variant is handed
/// out, so neither the variants nor the numbers depend on how many threads
/// there are. In wasm there are none and the same lines are parsed one
/// after another; a batch there holds one line.
#[cfg(not(target_family = "wasm"))]
fn parse_lines(lines: &mut [BatchLine], rules: &ParseRules<'_>) {
    use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

    lines.par_iter_mut().for_each(|line| line.parse(rules));
}

/// The lines of a batch, each parsed into the variant of its own line, one
/// after another, which is what wasm does: it has no threads.
#[cfg(target_family = "wasm")]
fn parse_lines(lines: &mut [BatchLine], rules: &ParseRules<'_>) {
    for line in lines {
        line.parse(rules);
    }
}

/// A reader over a VCF, which gives its variants one at a time.
///
/// It holds the individuals of the file and the names of the chromosomes
/// it has given so far, each with its number.
pub struct VcfReader<R: BufRead + Send> {
    source: VcfSource<R>,
    options: VcfOptions,
    individuals: Vec<String>,
    chroms: ChromTable,
    needs: Needs,
    /// The lines that were read together and parsed together, each with its
    /// variant. The ones after `filled` are the lines of the batches
    /// before, kept for their buffers.
    batch: Vec<BatchLine>,
    /// How many lines of the batch were read from the source.
    filled: usize,
    /// Which of them is the next to be handed out.
    next: usize,
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
    /// The error of the line that could not be read, which is given after
    /// the lines of the batch that were read before it.
    line_error: Option<Error>,
    /// Whether the source has been read to its end or could not be read.
    source_done: bool,
    /// Whether the reader gave its last variant or an error. After either,
    /// every read gives no variant.
    finished: bool,
    /// Whether a parse of a batch was begun and did not come back, which is
    /// what a panic inside the parse leaves behind: the lines of that batch
    /// hold what the batch before them left in them, so the reader cannot
    /// go on.
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
    /// When the ploidy of the options is 0 or above [`MAX_PLOIDY`], when
    /// the source is not a VCF,
    /// when its header has not the nine first columns of a VCF with
    /// genotypes or no individual after them, when two individuals have the
    /// same name, and when the source cannot be read.
    pub fn new(source: R, options: VcfOptions) -> Result<VcfReader<R>> {
        if options.ploidy == 0 || options.ploidy > MAX_PLOIDY {
            return Err(Error::VcfPloidyOutOfRange {
                ploidy: options.ploidy,
                largest: MAX_PLOIDY,
            });
        }
        let source = WithFirstBytes::new(source, BYTES_LOOKED_AT)?;
        let gzipped = source.first.starts_with(&GZIP_BYTES);
        let mut source = if gzipped {
            VcfSource::Gzipped(BufReader::new(MultiGzDecoder::new(source)))
        } else {
            VcfSource::Plain(source)
        };
        let first_bytes = source.first_bytes()?;
        if !first_bytes.starts_with(b"#") {
            return Err(Error::NotAVcf {
                found: as_text(first_bytes),
            });
        }
        let mut reader = VcfReader {
            source,
            options,
            individuals: Vec::new(),
            chroms: ChromTable::new(),
            needs: Needs::ALL,
            batch: Vec::new(),
            filled: 0,
            next: 0,
            lines_per_batch: LINES_PER_BATCH,
            bytes_per_batch: BYTES_PER_BATCH,
            #[cfg(test)]
            batches_filled: 0,
            line_number: 0,
            line_error: None,
            source_done: false,
            finished: false,
            parsing: false,
            #[cfg(test)]
            panic_at_line: None,
        };
        reader.read_header()?;
        Ok(reader)
    }

    /// It skips the `##` lines and takes the individuals from the `#CHROM`
    /// line, which it leaves consumed, so that the next line read is the
    /// first data line.
    ///
    /// The line it reads into is its own: the header is read once, when the
    /// reader is built, and the lines of the batch are not there yet.
    fn read_header(&mut self) -> Result<()> {
        let mut line = String::new();
        loop {
            line.clear();
            let number = next_line_number(self.line_number);
            let read = self.source.read_line(&mut line).map_err(|error| {
                if is_not_text(&error) {
                    Error::VcfHeader {
                        problem: format!(
                            "the bytes of its line {number} are not valid UTF-8, and a \
                             VCF is text"
                        ),
                    }
                } else {
                    Error::Io(error)
                }
            })?;
            if read == 0 {
                return Err(Error::VcfHeader {
                    problem: "it has no #CHROM line".to_string(),
                });
            }
            self.line_number = number;
            let text = without_the_line_end(&line);
            if text.starts_with("##") {
                continue;
            }
            self.individuals = individuals_of(text, self.line_number)?;
            return Ok(());
        }
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
        VcfReader::new(BufReader::new(file), options)
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

/// Whether the error of an input is bytes that are not valid UTF-8, which
/// is what `read_line` gives for a line of a VCF that is not text.
fn is_not_text(error: &std::io::Error) -> bool {
    error.kind() == std::io::ErrorKind::InvalidData
}

/// The line without the `\n` or the `\r\n` it ends in. The genotype of the
/// last individual is the one that would carry the `\r`.
fn without_the_line_end(line: &str) -> &str {
    let line = line.strip_suffix('\n').unwrap_or(line);
    line.strip_suffix('\r').unwrap_or(line)
}

/// The number of the next line. A file of `u64::MAX` lines cannot be
/// written, so the saturation is not reached and no error is worth adding
/// for it.
fn next_line_number(line_number: u64) -> u64 {
    line_number.saturating_add(1)
}

/// The first bytes of a source as text, for the message that says it is not
/// a VCF: the ones that are not printable go in as their number.
fn as_text(bytes: &[u8]) -> String {
    if bytes.is_empty() {
        return "nothing: it holds no byte".to_string();
    }
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

/// Whether the FILTER column says that the variant passed: `PASS`, or a
/// dot, which says that no filter was applied to it.
fn passed(filter: &str) -> bool {
    filter == "PASS" || filter == MISSING_VALUE
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
fn parse_quality(text: &str, line: u64) -> Result<Option<f32>> {
    if text == MISSING_VALUE {
        return Ok(None);
    }
    text.parse().map(Some).map_err(|_| Error::VcfDataLine {
        line,
        place: VcfPlace::Column("QUAL"),
        problem: format!("`{text}` is not a quality"),
    })
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

/// The alleles of the variant, written over the strings of `spare`, the
/// ones the reader took out of the variant before this read.
///
/// A variant with more alleles than the one before it would allocate a
/// string if the strings that are left over were dropped instead of kept:
/// a second pass over `many.vcf` with every field asked for allocated 54
/// times, once for each of its 54 variants with two alternative alleles.
/// Section 1 of `docs/architecture.md` asks for none.
fn fill_alleles(
    alleles: &mut Vec<String>,
    spare: &mut Vec<String>,
    reference: &str,
    alternatives: &str,
) {
    for text in allele_texts(reference, alternatives) {
        let mut allele = spare.pop().unwrap_or_default();
        allele.clear();
        allele.push_str(text);
        alleles.push(allele);
    }
}

/// One allele of a genotype: a number of the alleles the variant declares,
/// or [`MISSING_ALLELE`] for a dot.
fn parse_allele(text: &str, num_alleles: usize, line: u64, individual: &str) -> Result<i8> {
    if text == MISSING_VALUE {
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
    for byte in text.as_bytes() {
        let Some(digit) = char::from(*byte).to_digit(10) else {
            return Err(wrong(format!(
                "`{text}` is not an allele number, which is a run of digits"
            )));
        };
        // Once the number is above the largest allele the answer is the
        // same whatever its other digits are, so it stops growing there
        // and the two operations cannot overflow.
        number = number.saturating_mul(10).saturating_add(digit);
    }
    let Ok(allele) = i8::try_from(number) else {
        return Err(wrong(format!(
            "the allele `{text}` is above {MAX_ALLELE}, the largest allele popnei holds"
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

/// The alleles of the genotype of one individual, appended to `gts`. A
/// genotype written as a single dot is a missing genotype of the ploidy of
/// the file, and any other number of alleles than the ploidy is an error.
fn fill_genotype(
    gts: &mut Vec<i8>,
    text: &str,
    ploidy: usize,
    num_alleles: usize,
    line: u64,
    individual: &str,
) -> Result<()> {
    // VCF 4.4 lets a genotype start with its separator, `/0/1`.
    let text = text.strip_prefix(['/', '|']).unwrap_or(text);
    if text == MISSING_VALUE {
        gts.extend(std::iter::repeat_n(MISSING_ALLELE, ploidy));
        return Ok(());
    }
    // The alleles are read in one pass over the text, and what was pushed
    // is taken off again when the genotype turns out to be of another
    // ploidy or one of its alleles is wrong, so that a genotype leaves
    // either all of its alleles in `gts` or none.
    let start = gts.len();
    for allele in text.split(['/', '|']) {
        match parse_allele(allele, num_alleles, line, individual) {
            Ok(allele) => gts.push(allele),
            Err(error) => {
                gts.truncate(start);
                return Err(error);
            }
        }
    }
    // `gts` only grew in the loop above, so the subtraction does not
    // saturate: it is how many alleles this genotype pushed.
    let found = gts.len().saturating_sub(start);
    if found != ploidy {
        gts.truncate(start);
        return Err(Error::VcfGenotypePloidy {
            line,
            individual: individual.to_string(),
            found,
            expected: ploidy,
        });
    }
    Ok(())
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

/// The genotypes of every individual of the line, appended to `gts`: the
/// value of the key `GT` of each column, which an individual that drops its
/// last values still has.
fn fill_genotypes<'a>(
    gts: &mut Vec<i8>,
    columns: &mut impl Iterator<Item = &'a str>,
    gt_index: usize,
    individuals: &[String],
    ploidy: usize,
    num_alleles: usize,
    line: u64,
) -> Result<()> {
    for (read_so_far, individual) in individuals.iter().enumerate() {
        let Some(column) = columns.next() else {
            return Err(Error::VcfDataLine {
                line,
                place: VcfPlace::Line,
                problem: format!(
                    "it has the columns of {read_so_far} individuals and the header has {count}",
                    count = individuals.len(),
                ),
            });
        };
        let Some(genotype) = column.split(':').nth(gt_index) else {
            return Err(Error::VcfDataLine {
                line,
                place: VcfPlace::Individual(individual.clone()),
                problem: format!("`{column}` has no value where the FORMAT has GT"),
            });
        };
        fill_genotype(gts, genotype, ploidy, num_alleles, line, individual)?;
    }
    let left_over = columns.count();
    if left_over != 0 {
        return Err(Error::VcfDataLine {
            line,
            place: VcfPlace::Line,
            problem: format!(
                "it has {left_over} columns more than the {count} individuals of the header",
                count = individuals.len(),
            ),
        });
    }
    Ok(())
}

/// What the parse of a data line needs to know, which is the same for
/// every line of a file: the reader hands one to each line it parses.
struct ParseRules<'a> {
    options: VcfOptions,
    needs: Needs,
    individuals: &'a [String],
    /// The line whose parse panics. No VCF makes the parse panic, and this
    /// is how the test of what a reader does after a panic in its parse
    /// makes one happen; nothing outside the tests can set it.
    #[cfg(test)]
    panic_at_line: Option<u64>,
}

/// The variant of the data line `text`, the line `number` of the file,
/// into `var`, and false when the line gives no variant: an empty line, or
/// one whose FILTER failed when the options leave those out.
///
/// The name of the chromosome goes into `chrom_name` and its number is not
/// given here: the reader gives it when it hands the variant out, so that
/// the numbers follow the order of the file and not the order in which the
/// lines were parsed. The alleles are written over the strings of
/// `spare_alleles`, the ones of the variants read before.
///
/// `var` is cleared by the caller, which is also what puts the strings of
/// the alleles it held into `spare_alleles`.
fn parse_data_line(
    text: &str,
    number: u64,
    rules: &ParseRules<'_>,
    var: &mut Variant,
    chrom_name: &mut String,
    spare_alleles: &mut Vec<String>,
) -> Result<bool> {
    let ParseRules {
        options,
        needs,
        individuals,
        ..
    } = rules;
    if text.is_empty() {
        return Ok(false);
    }
    // The seven columns up to the FILTER are taken as text and read only
    // when the variant is given: a line that is skipped costs no parsing,
    // and the name of its chromosome gets no number.
    let mut columns = text.split('\t');
    let chrom_text = next_column(&mut columns, "CHROM", number)?;
    let pos_text = next_column(&mut columns, "POS", number)?;
    let id_text = next_column(&mut columns, "ID", number)?;
    let reference_text = next_column(&mut columns, "REF", number)?;
    let alternatives_text = next_column(&mut columns, "ALT", number)?;
    let quality_text = next_column(&mut columns, "QUAL", number)?;
    let filter_text = next_column(&mut columns, "FILTER", number)?;
    if options.only_passed && !passed(filter_text) {
        return Ok(false);
    }

    var.pos = parse_position(pos_text, number)?;
    chrom_name.push_str(chrom_text);
    var.filled = Needs::CHROM_POS;
    if needs.contains(Needs::ID) {
        if id_text != MISSING_VALUE {
            var.id.push_str(id_text);
        }
        var.filled |= Needs::ID;
    }
    // The alleles are counted for every variant that is given, to check the
    // allele numbers of its genotypes, also when the texts of the alleles
    // are not kept.
    let num_alleles = count_alleles(reference_text, alternatives_text, number)?;
    if needs.contains(Needs::ALLELES) {
        fill_alleles(
            &mut var.alleles,
            spare_alleles,
            reference_text,
            alternatives_text,
        );
        var.filled |= Needs::ALLELES;
    }
    if needs.contains(Needs::QUAL) {
        var.qual = parse_quality(quality_text, number)?;
        var.filled |= Needs::QUAL;
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
    let first_individual = next_column(&mut columns, "individual", number)?;
    if needs.contains(Needs::GTS) {
        let mut columns = std::iter::once(first_individual).chain(columns);
        fill_genotypes(
            &mut var.gts,
            &mut columns,
            gt_index,
            individuals,
            options.ploidy,
            num_alleles,
            number,
        )?;
        var.filled |= Needs::GTS;
    }
    Ok(true)
}

/// What the reader gives when a line of the source could not be read: the
/// bytes that are not text are an error of that line, with its number,
/// since a VCF is text and the number of the line is what a user needs, and
/// anything else is an error of the input.
fn error_reading_a_line(error: std::io::Error, number: u64) -> Error {
    if is_not_text(&error) {
        Error::VcfDataLine {
            line: number,
            place: VcfPlace::Line,
            problem: "its bytes are not valid UTF-8, and a VCF is text".to_string(),
        }
    } else {
        Error::Io(error)
    }
}

impl<R: BufRead + Send> VcfReader<R> {
    /// How many lines the reader takes from the source before it parses
    /// them, which is [`LINES_PER_BATCH`] until this is called. The tests
    /// lower it, so that a file of a few hundred lines is read in several
    /// batches and so that one line at a time, which is what wasm reads, is
    /// read here too. A batch holds one line at least.
    #[cfg(test)]
    fn set_lines_per_batch(&mut self, lines: usize) {
        self.lines_per_batch = lines.max(1);
    }

    /// How many bytes of text the reader takes from the source before it
    /// parses what it read, which is [`BYTES_PER_BATCH`] until this is
    /// called. A batch holds one line at least, whatever this says.
    #[cfg(test)]
    fn set_bytes_per_batch(&mut self, bytes: usize) {
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

    /// It reads the next lines of the source into the batch, over the text
    /// and the buffers of the batch before, and parses them.
    ///
    /// A line that cannot be read ends the batch and its error is kept
    /// apart, to be given after the lines that were read before it.
    fn fill_batch(&mut self) {
        let VcfReader {
            source,
            options,
            individuals,
            needs,
            batch,
            filled,
            next,
            lines_per_batch,
            bytes_per_batch,
            line_number,
            line_error,
            source_done,
            parsing,
            #[cfg(test)]
            batches_filled,
            #[cfg(test)]
            panic_at_line,
            ..
        } = self;
        *next = 0;
        *filled = 0;
        #[cfg(test)]
        {
            *batches_filled = batches_filled.saturating_add(1);
        }
        // The text of the lines that were read, which bounds the batch
        // beside their number: one line of a file of many individuals is
        // where the memory of a reader would otherwise grow without a
        // bound.
        let mut bytes: usize = 0;
        while *filled < *lines_per_batch && bytes < *bytes_per_batch {
            if batch.len() <= *filled {
                batch.push(BatchLine::new());
            }
            let Some(line) = batch.get_mut(*filled) else {
                break;
            };
            line.text.clear();
            let number = next_line_number(*line_number);
            match source.read_line(&mut line.text) {
                Ok(0) => {
                    *source_done = true;
                    break;
                }
                Ok(read) => {
                    *line_number = number;
                    line.number = number;
                    // A batch holds `lines_per_batch` lines at most, so the
                    // count does not reach the largest `usize`, and the
                    // bytes of a batch stop growing at the line that
                    // reaches the bound.
                    *filled = filled.saturating_add(1);
                    bytes = bytes.saturating_add(read);
                }
                Err(error) => {
                    *line_number = number;
                    *line_error = Some(error_reading_a_line(error, number));
                    *source_done = true;
                    break;
                }
            }
        }
        let rules = ParseRules {
            options: *options,
            needs: *needs,
            individuals,
            #[cfg(test)]
            panic_at_line: *panic_at_line,
        };
        if let Some(lines) = batch.get_mut(..*filled) {
            // A panic of a worker unwinds through here, and what says so
            // afterwards is this flag, which is set until the parse comes
            // back.
            *parsing = true;
            parse_lines(lines, &rules);
            *parsing = false;
        }
    }

    /// The next variant of the file, from the batch that was parsed, and
    /// from the next batch when that one is spent: the empty lines and,
    /// when the options ask for it, the variants that failed a filter are
    /// the lines of the batch that give no variant.
    fn next_variant(&mut self, var: &mut Variant) -> Result<bool> {
        if self.parsing {
            return Err(Error::VcfParseNotFinished {
                line: self.line_number,
            });
        }
        if self.finished {
            return Ok(false);
        }
        loop {
            while self.next < self.filled {
                let VcfReader {
                    chroms,
                    batch,
                    next,
                    ..
                } = self;
                let index = *next;
                // A batch holds `lines_per_batch` lines at most, so the
                // count does not reach the largest `usize`.
                *next = index.saturating_add(1);
                let Some(line) = batch.get_mut(index) else {
                    break;
                };
                // What the line gave is taken out of it: the line is left
                // with the variant of the consumer, which the next batch
                // written into this place clears, and with nothing to give
                // again.
                match std::mem::replace(&mut line.outcome, LineOutcome::NoVariant) {
                    LineOutcome::NoVariant => {}
                    LineOutcome::Variant => {
                        // The name of the chromosome gets its number here,
                        // when the variant is handed out, so that the
                        // numbers follow the order of the file whatever the
                        // threads did.
                        let chrom = chroms.intern(&line.chrom_name);
                        line.var.chrom = chrom;
                        std::mem::swap(var, &mut line.var);
                        return Ok(true);
                    }
                    LineOutcome::Wrong(error) => return Err(error),
                }
            }
            if let Some(error) = self.line_error.take() {
                return Err(error);
            }
            if self.source_done {
                self.finished = true;
                return Ok(false);
            }
            self.fill_batch();
        }
    }
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
            .field("chroms", &self.chroms)
            .field("line_number", &self.line_number)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl<R: BufRead + Send> VariantReader for VcfReader<R> {
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
        // The variant of the consumer is swapped with the one of the line
        // that is handed out, so what it held goes back into the batch and
        // is written over there, the strings of its alleles among them.
        match self.next_variant(var) {
            Ok(true) => Ok(true),
            Ok(false) => {
                var.clear();
                Ok(false)
            }
            Err(error) => {
                self.finished = true;
                var.clear();
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

    fn set_needs(&mut self, needs: Needs) {
        if needs == self.needs {
            return;
        }
        self.needs = needs;
        // The lines of the batch that were read and not handed out yet were
        // parsed for the fields that were asked for before, and the change
        // holds from the next read on, so they are parsed again: their text
        // is still in the batch and the parse writes over the same buffers.
        // A line that could not be read is not one of them, and its error
        // waits apart from the batch.
        let VcfReader {
            options,
            individuals,
            needs,
            batch,
            filled,
            next,
            #[cfg(test)]
            panic_at_line,
            ..
        } = self;
        let rules = ParseRules {
            options: *options,
            needs: *needs,
            individuals,
            #[cfg(test)]
            panic_at_line: *panic_at_line,
        };
        if let Some(lines) = batch.get_mut(*next..*filled) {
            parse_lines(lines, &rules);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use super::{BYTES_PER_BATCH, MAX_PLOIDY, MISSING_VALUE, VcfOptions, VcfPlace, VcfReader};
    use crate::error::{Error, Result};
    use crate::variant::{MISSING_ALLELE, Needs, Variant, VariantReader};

    /// The panic that the test of a reader whose parse did not come back
    /// injects into the parse of one line. It is the only way to make the
    /// parse panic: no VCF does it, and every error of a line is a value
    /// that the line carries back.
    pub(super) fn panic_if_the_test_asked_for_it(number: u64, rules: &super::ParseRules<'_>) {
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

    #[test]
    fn the_individuals_of_a_plain_and_of_a_gzipped_vcf_are_read() {
        for name in ["cases.vcf", "cases.vcf.gz"] {
            let reader = VcfReader::from_path(&reference_vcf(name), VcfOptions::default())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
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
    fn a_gzipped_source_cut_in_the_middle_of_a_member_is_refused() {
        let bytes = std::fs::read(reference_vcf("cases.vcf.gz")).unwrap();
        // The first member of this file is its header and the second its
        // four variants, so the cut is inside the second one.
        let cut = bytes.len().saturating_sub(20);
        let bytes = bytes.get(..cut).unwrap_or_default().to_vec();
        let mut reader = VcfReader::new(Cursor::new(bytes), VcfOptions::default()).unwrap();
        assert_eq!(reader.individuals(), ["ind1", "ind2", "ind3"]);
        let error = rows_of(&mut reader).unwrap_err();
        assert!(
            matches!(error, Error::Io(_)),
            "a file cut in the middle of a gzip member gives {error}"
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
    fn a_read_after_an_error_and_after_the_last_variant_gives_no_variant() {
        let mut reader = reader_over(
            &vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]),
            VcfOptions::default(),
        );
        let mut var = Variant::new();
        assert!(reader.read_variant(&mut var).unwrap());
        assert!(!reader.read_variant(&mut var).unwrap());
        assert!(!reader.read_variant(&mut var).unwrap());

        let mut reader = reader_over(
            &vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1"]),
            VcfOptions::default(),
        );
        assert!(reader.read_variant(&mut var).is_err());
        assert!(!reader.read_variant(&mut var).unwrap());
        assert!(!reader.read_variant(&mut var).unwrap());
        assert_eq!(var.filled, Needs::empty());
    }

    /// A reader over bytes held in memory, which is how the tests give a
    /// VCF of their own.
    fn reader_over(vcf: &str, options: VcfOptions) -> VcfReader<Cursor<Vec<u8>>> {
        match VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), options) {
            Ok(reader) => reader,
            Err(error) => panic!("the reader was not built: {error}"),
        }
    }

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

    /// One variant as the tables of "How it is verified" of
    /// `docs/specs/io_vcf.md` give it, with the name of the chromosome in
    /// the place of the number the reader gave it.
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

    /// The variant a reader filled, with the name that its table of
    /// chromosomes gives to the number of that variant.
    fn row_of(var: &Variant, reader: &impl VariantReader) -> Row {
        Row {
            chrom: reader
                .chroms()
                .name(var.chrom)
                .unwrap_or("no name")
                .to_string(),
            pos: var.pos,
            id: var.id.clone(),
            alleles: var.alleles.clone(),
            qual: var.qual,
            gts: var.gts.clone(),
        }
    }

    /// Every variant a reader gives, until it has no more or it fails.
    fn rows_of(reader: &mut impl VariantReader) -> Result<Vec<Row>> {
        let mut var = Variant::new();
        let mut rows = Vec::new();
        while reader.read_variant(&mut var)? {
            rows.push(row_of(&var, reader));
        }
        Ok(rows)
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
        let mut reader = match VcfReader::from_path(&reference_vcf(name), options) {
            Ok(reader) => reader,
            Err(error) => panic!("{name}: {error}"),
        };
        match rows_of(&mut reader) {
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

    /// The options of a file of another ploidy, or of one read with every
    /// variant given.
    fn options(ploidy: usize, only_passed: bool) -> VcfOptions {
        VcfOptions {
            ploidy,
            only_passed,
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
        let mut reader = reader_over(&vcf, VcfOptions::default());
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!((var.pos, var.chrom), (20, 0));
        assert_eq!(var.filled, Needs::ALL);
        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!((var.pos, var.chrom), (30, 1));

        assert!(!reader.read_variant(&mut var).unwrap());
        assert_eq!(reader.chroms().len(), 2);
        assert_eq!(reader.chroms().name(0), Some("chr1"));
        assert_eq!(reader.chroms().name(1), Some("chr9"));
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
    fn a_haploid_genotype_with_a_ploidy_of_two_is_refused_after_the_variants_before_it() {
        let vcf = vcf_of(&[
            "chr1 100 . A T . PASS . GT 0/0 0/1 1/1",
            "chr1 200 . A T . PASS . GT 1 0 1",
        ]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.pos, 100);

        let error = reader.read_variant(&mut var).unwrap_err();
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
        reader.set_needs(Needs::GTS);
        let mut var = Variant::new();

        let error = reader.read_variant(&mut var).unwrap_err();
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
    fn a_vcf_with_a_header_and_no_variant_gives_no_variant_and_no_error() {
        let mut reader = reader_over(HEADER, VcfOptions::default());
        let mut var = Variant::new();
        assert!(!reader.read_variant(&mut var).unwrap());
        assert!(!reader.read_variant(&mut var).unwrap());
        assert_eq!(var.filled, Needs::empty());
        assert!(var.gts.is_empty());
        assert!(reader.chroms().is_empty());
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
    fn with_the_genotypes_alone_the_id_and_the_alleles_are_not_filled() {
        let vcf = vcf_of(&["chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::GTS);
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.filled, Needs::GTS | Needs::CHROM_POS);
        assert_eq!(var.gts, [0, 0, 0, 1, 1, 1]);
        assert_eq!(var.pos, 100);
        assert_eq!(reader.chroms().name(var.chrom), Some("chr1"));
        assert!(var.alleles.is_empty());
        assert!(var.id.is_empty());
        assert_eq!(var.qual, None);
    }

    #[test]
    fn with_the_id_and_the_alleles_alone_the_genotypes_are_not_filled() {
        let vcf = vcf_of(&["chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1"]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        reader.set_needs(Needs::ID | Needs::ALLELES);
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.filled, Needs::ID | Needs::ALLELES | Needs::CHROM_POS);
        assert_eq!(var.id, "rs1");
        assert_eq!(var.alleles, ["A", "T"]);
        assert_eq!(var.pos, 100);
        assert!(var.gts.is_empty());
        assert_eq!(var.qual, None);
    }

    /// The variants of a VCF read with the id and the alleles asked for
    /// and not the genotypes, or the error the reader stops at.
    fn read_without_the_genotypes(vcf: &str) -> Result<Vec<Row>> {
        let mut reader = reader_over(vcf, VcfOptions::default());
        reader.set_needs(Needs::ID | Needs::ALLELES);
        rows_of(&mut reader)
    }

    #[test]
    fn a_data_line_whose_bytes_are_not_text_is_refused_with_its_number() {
        let mut vcf = vcf_of(&["chr1 100 . A T . PASS . GT 0/0 0/1 1/1"]).into_bytes();
        vcf.extend_from_slice(b"chr1\t200\t.\t\xffA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\n");
        let mut reader = VcfReader::new(Cursor::new(vcf), VcfOptions::default()).unwrap();
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.pos, 100);

        let error = reader.read_variant(&mut var).unwrap_err();
        let Error::VcfDataLine {
            line,
            place,
            problem,
        } = error
        else {
            panic!("the error is {error}");
        };
        assert_eq!((line, place), (5, VcfPlace::Line));
        assert!(problem.contains("UTF-8"), "{problem}");
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
            vec![
                row("chr1", 100, "rs1", &["A", "T"], None, &[]),
                row("chr1", 200, "rs2", &["A", "T"], None, &[]),
            ]
        );
    }

    #[test]
    fn a_reader_asked_for_the_genotypes_alone_empties_the_alleles_it_filled_before() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 A T 29.5 PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T 29.5 PASS . GT 0/0 0/1 1/1",
        ]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.alleles, ["A", "T"]);
        assert_eq!(var.id, "rs1");
        assert_eq!(var.filled, Needs::ALL);

        reader.set_needs(Needs::GTS);
        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.filled, Needs::GTS | Needs::CHROM_POS);
        assert!(var.alleles.is_empty());
        assert!(var.id.is_empty());
        assert_eq!(var.qual, None);
    }

    #[test]
    fn the_buffers_of_a_variant_are_written_over_from_one_variant_to_the_next() {
        let vcf = vcf_of(&[
            "chr1 100 rs1 AAAA TTTT . PASS . GT 0/0 0/1 1/1",
            "chr1 200 rs2 A T,G . PASS . GT 0/0 0/1 1/2",
            "chr1 300 . A T . PASS . GT 0/0 0/1 1/1",
        ]);
        let mut reader = reader_over(&vcf, VcfOptions::default());
        let mut var = Variant::new();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.alleles, ["AAAA", "TTTT"]);
        let gts = var.gts.capacity();

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.alleles, ["A", "T", "G"]);
        assert_eq!(var.id, "rs2");

        assert!(reader.read_variant(&mut var).unwrap());
        assert_eq!(var.alleles, ["A", "T"]);
        assert_eq!(var.id, "");
        assert_eq!(var.gts.capacity(), gts);
    }

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
    /// its chromosome, read inside a rayon pool of `threads` threads and in
    /// batches of `lines_per_batch` lines.
    ///
    /// The pool is built here and is not rayon's global one, which has one
    /// thread per core of the machine: `install` runs the reader on this
    /// one instead. rayon is a dependency of the targets that are not wasm,
    /// so this and what uses it are compiled for those alone.
    #[cfg(not(target_family = "wasm"))]
    fn many_vcf_read_in_a_pool(threads: usize, lines_per_batch: usize) -> Vec<(u32, Row)> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();
        pool.install(|| {
            let mut reader =
                match VcfReader::from_path(&reference_vcf("many.vcf"), options(MANY_PLOIDY, false))
                {
                    Ok(reader) => reader,
                    Err(error) => panic!("many.vcf: {error}"),
                };
            reader.set_lines_per_batch(lines_per_batch);
            let mut var = Variant::new();
            let mut variants = Vec::new();
            while match reader.read_variant(&mut var) {
                Ok(read) => read,
                Err(error) => panic!("many.vcf: the reader stopped at {error}"),
            } {
                variants.push((var.chrom, row_of(&var, &reader)));
            }
            variants
        })
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn many_vcf_gives_the_same_variants_in_a_pool_of_one_thread_and_in_one_of_four() {
        // 64 lines a batch, and the file has 500 data lines, so the pool of
        // four threads parses eight batches and every thread parses lines
        // of each.
        let one_thread = many_vcf_read_in_a_pool(1, 64);
        let four_threads = many_vcf_read_in_a_pool(4, 64);
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
    fn a_reader_whose_parse_did_not_come_back_gives_an_error_and_no_variant() {
        // Eight data lines, in batches of four, and the parse of the third
        // line of the second batch panics: the four variants of the first
        // batch are given, and then the read that fills the second batch
        // panics inside the parse and unwinds through the reader.
        let lines: Vec<String> = (1..=8)
            .map(|pos| format!("chr1 {pos}00 . A T . PASS . GT 0/0 0/1 1/1"))
            .collect();
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let mut reader = reader_over(&vcf_of(&lines), VcfOptions::default());
        reader.set_lines_per_batch(4);
        // The three lines of the header come first, so the sixth data line
        // is the line 9 of the file.
        reader.panic_at_line(9);
        let mut var = Variant::new();

        let mut read = 0;
        let hook = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let went_on = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            while reader.read_variant(&mut var).unwrap() {
                read += 1;
            }
        }));
        std::panic::set_hook(hook);
        assert!(went_on.is_err(), "the parse did not panic");
        assert_eq!(read, 4);

        // The lines of the batch that was being parsed were never parsed,
        // and a reader that went on would give the variants of the ones
        // that were and drop the others without a word.
        let error = reader.read_variant(&mut var).unwrap_err();
        let Error::VcfParseNotFinished { line } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(line, 11);
        // Every read after it gives the same error and not the false of a
        // file that was read to its end, which a consumer would take for a
        // VCF that ends there.
        let again = reader.read_variant(&mut var).unwrap_err();
        assert!(
            matches!(again, Error::VcfParseNotFinished { line: 11 }),
            "the second error is {again}"
        );
        assert_eq!(var.filled, Needs::empty());
    }

    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn a_wrong_line_of_a_later_batch_comes_after_the_variants_that_were_read_before_it() {
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
            let mut reader = reader_over(&vcf, VcfOptions::default());
            reader.set_lines_per_batch(64);
            let mut var = Variant::new();
            for pos in 1..=200 {
                assert!(reader.read_variant(&mut var).unwrap(), "the variant {pos}");
                assert_eq!(var.pos, pos);
            }

            let error = reader.read_variant(&mut var).unwrap_err();
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
            assert!(!reader.read_variant(&mut var).unwrap());
        });
    }

    #[test]
    fn a_bound_of_bytes_smaller_than_the_file_cuts_the_batches_and_changes_no_result() {
        let expected = rows_of_file("many.vcf", VcfOptions::default());
        assert_eq!(expected.len(), 475);
        // `many.vcf` has 500 data lines of about 230 bytes each. A bound of
        // one byte holds one line in every batch, which is the line a batch
        // holds whatever the bound says; the batch of the last line stops
        // at the bound without seeing the end of the file, so one more
        // batch is filled, which reads no line: 501. The bound of the code,
        // 8 MiB, takes the 500 lines and the end of the file in one.
        for (bytes, batches) in [(1, 501), (BYTES_PER_BATCH, 1)] {
            let mut reader =
                VcfReader::from_path(&reference_vcf("many.vcf"), VcfOptions::default())
                    .unwrap_or_else(|error| panic!("many.vcf: {error}"));
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
        let mut reader = VcfReader::from_path(&reference_vcf("many.vcf"), VcfOptions::default())
            .unwrap_or_else(|error| panic!("many.vcf: {error}"));
        reader.set_lines_per_batch(1);
        let one_line_at_a_time = rows_of(&mut reader).unwrap();
        assert_eq!(
            one_line_at_a_time,
            rows_of_file("many.vcf", VcfOptions::default())
        );
        assert_eq!(one_line_at_a_time.len(), 475);
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
}
