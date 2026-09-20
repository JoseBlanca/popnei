//! The VCF reader: it gives the variants of a VCF, plain or gzipped, one at
//! a time.
//!
//! [`VcfReader::new`] takes any source of bytes and reads the header of the
//! VCF, so the individuals are known before a variant is, and
//! [`VcfReader::from_path`] does the same for a caller that has a path. The
//! reader is a [`VariantReader`]: the consumer owns one
//! [`Variant`](crate::variant::Variant) and lends it to `read_variant` again
//! and again.
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
use crate::variant::{ChromTable, Needs, Variant, VariantReader};

/// The ploidy a VCF is read with when the caller asks for no other, the
/// ploidy of a diploid organism. pyNei has no such argument and reports 2
/// for every file it reads.
pub const DEFAULT_PLOIDY: usize = 2;

/// Whether a VCF is read with the variants that failed a filter left out,
/// when the caller asks for no other. The owner decided in September 2026
/// that the FILTER column is honoured and that this is the default; pyNei
/// ignores that column and gives every variant.
pub const DEFAULT_ONLY_PASSED: bool = true;

/// The two bytes every gzipped file starts with.
const GZIP_BYTES: [u8; 2] = [0x1f, 0x8b];

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
    Column(String),
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

/// The bytes of the VCF, decompressed when the source was gzipped.
enum VcfSource<R: BufRead> {
    /// The source as it was given.
    Plain(R),
    /// The source through the decoder of gzip, which goes on to the next
    /// member of the file when one ends.
    Gzipped(BufReader<MultiGzDecoder<R>>),
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
    fn read_line(&mut self, line: &mut String) -> std::io::Result<usize> {
        match self {
            VcfSource::Plain(source) => source.read_line(line),
            VcfSource::Gzipped(source) => source.read_line(line),
        }
    }
}

/// A reader over a VCF, which gives its variants one at a time.
///
/// It holds the individuals and the names of the chromosomes it has seen,
/// and the buffer of one line, which is refilled for every line of the
/// file.
pub struct VcfReader<R: BufRead + Send> {
    source: VcfSource<R>,
    options: VcfOptions,
    individuals: Vec<String>,
    chroms: ChromTable,
    needs: Needs,
    /// The line being read, without its end of line.
    line: String,
    /// The number of that line in the file, counted from 1 with the lines
    /// of the header.
    line_number: u64,
    /// Whether the source has been read to its end or gave an error. After
    /// either, every read gives no variant.
    finished: bool,
}

impl<R: BufRead + Send> VcfReader<R> {
    /// The reader over `source`, the VCF, gzipped or not, whose header it
    /// reads: the individuals are known when it returns.
    ///
    /// # Errors
    ///
    /// When the ploidy of the options is 0, when the source is not a VCF,
    /// when its header has not the nine first columns of a VCF with
    /// genotypes or no individual after them, when two individuals have the
    /// same name, and when the source cannot be read.
    pub fn new(source: R, options: VcfOptions) -> Result<VcfReader<R>> {
        if options.ploidy == 0 {
            return Err(Error::VcfPloidyIsZero);
        }
        let mut source = source;
        let gzipped = source.fill_buf()?.starts_with(&GZIP_BYTES);
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
            line: String::new(),
            line_number: 0,
            finished: false,
        };
        reader.read_header()?;
        Ok(reader)
    }

    /// It skips the `##` lines and takes the individuals from the `#CHROM`
    /// line, which it leaves consumed, so that the next line read is the
    /// first data line.
    fn read_header(&mut self) -> Result<()> {
        loop {
            self.line.clear();
            let read = self.source.read_line(&mut self.line)?;
            if read == 0 {
                return Err(Error::VcfHeader {
                    problem: "it has no #CHROM line".to_string(),
                });
            }
            self.line_number = next_line_number(self.line_number);
            let text = without_the_line_end(&self.line);
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
    /// When the file cannot be opened or read, and everything
    /// [`VcfReader::new`] fails with.
    pub fn from_path(path: &Path, options: VcfOptions) -> Result<Self> {
        let file = File::open(path)?;
        VcfReader::new(BufReader::new(file), options)
    }
}

/// The individuals of the `#CHROM` line, whose nine first columns have to
/// be the nine of a VCF with genotypes.
fn individuals_of(chrom_line: &str, line_number: u64) -> Result<Vec<String>> {
    let columns: Vec<&str> = chrom_line.split('\t').collect();
    for (found, expected) in columns.iter().zip(FIRST_COLUMNS) {
        if *found != expected {
            return Err(Error::VcfHeader {
                problem: format!(
                    "its line {line_number} has `{found}` where a VCF with genotypes \
                     has `{expected}`; the nine first columns are {first}",
                    first = FIRST_COLUMNS.join(" "),
                ),
            });
        }
    }
    let Some(individuals) = columns.get(FIRST_COLUMNS.len()..) else {
        return Err(Error::VcfHeader {
            problem: format!(
                "its line {line_number} has {count} columns, and a VCF with genotypes \
                 has the nine {first} and one column per individual after them",
                count = columns.len(),
                first = FIRST_COLUMNS.join(" "),
            ),
        });
    };
    if individuals.is_empty() {
        return Err(Error::VcfHeader {
            problem: format!(
                "its line {line_number} has a FORMAT column and no individual after it"
            ),
        });
    }
    let mut seen = HashSet::with_capacity(individuals.len());
    for name in individuals {
        if !seen.insert(*name) {
            return Err(Error::VcfHeader {
                problem: format!("two individuals of its line {line_number} are called `{name}`"),
            });
        }
    }
    Ok(individuals.iter().map(|name| (*name).to_string()).collect())
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
    let mut text = String::new();
    for byte in bytes.iter().take(16) {
        if byte.is_ascii_graphic() || *byte == b' ' {
            text.push(char::from(*byte));
        } else {
            text.push_str(&format!("\\x{byte:02x}"));
        }
    }
    format!("`{text}`")
}

impl<R: BufRead + Send> VariantReader for VcfReader<R> {
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
        // The data lines are the next commit, task 3.2 of the plan. Until
        // then the reader reads the header and gives no variant.
        var.clear();
        if self.finished {
            return Ok(false);
        }
        self.finished = true;
        Ok(false)
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
        self.needs = needs;
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::{VcfOptions, VcfReader};
    use crate::error::Error;
    use crate::variant::VariantReader;

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
        assert!(problem.contains("#CHROM"), "{problem}");
    }

    #[test]
    fn a_ploidy_of_zero_is_refused() {
        let options = VcfOptions {
            ploidy: 0,
            only_passed: true,
        };
        let error = error_of(HEADER, options);
        assert!(
            matches!(error, Error::VcfPloidyIsZero),
            "the error is {error}"
        );
    }
}
