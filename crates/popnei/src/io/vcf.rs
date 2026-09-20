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

/// The two bytes every gzipped file starts with.
const GZIP_BYTES: [u8; 2] = [0x1f, 0x8b];

/// What a column of a VCF holds when it has no value: the id of a variant
/// with no id, an ALT with no alternative allele, an allele that was not
/// called.
const MISSING_VALUE: &str = ".";

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
    if bytes.is_empty() {
        return "nothing: it holds no byte".to_string();
    }
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
        place: VcfPlace::Column("POS".to_string()),
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
        place: VcfPlace::Column("QUAL".to_string()),
        problem: format!("`{text}` is not a quality"),
    })
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

/// The alleles of the variant, written over the strings that `alleles`
/// holds already: once they have grown, a million variants cost no
/// allocation, which is what section 1 of `docs/architecture.md` asks of a
/// reader. `num_alleles` is how many texts there are.
fn fill_alleles(
    alleles: &mut Vec<String>,
    reference: &str,
    alternatives: &str,
    num_alleles: usize,
) {
    for (position, text) in allele_texts(reference, alternatives).enumerate() {
        match alleles.get_mut(position) {
            Some(allele) => {
                allele.clear();
                allele.push_str(text);
            }
            None => alleles.push(text.to_string()),
        }
    }
    alleles.truncate(num_alleles);
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
    let Ok(number) = text.parse::<u32>() else {
        return Err(wrong(format!("`{text}` is not an allele number")));
    };
    let Ok(allele) = i8::try_from(number) else {
        return Err(wrong(format!(
            "the allele {number} is above {MAX_ALLELE}, the largest allele popnei holds"
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
    // The alleles are counted before any of them is read, so that a
    // genotype of another ploidy leaves nothing in `gts`.
    let found = text.split(['/', '|']).count();
    if found != ploidy {
        return Err(Error::VcfGenotypePloidy {
            line,
            individual: individual.to_string(),
            found,
            expected: ploidy,
        });
    }
    for allele in text.split(['/', '|']) {
        gts.push(parse_allele(allele, num_alleles, line, individual)?);
    }
    Ok(())
}

/// The genotypes of every individual of the line, appended to `gts`: the
/// value of the key `GT` of each column, which an individual that drops its
/// last values still has.
fn fill_genotypes<'a>(
    gts: &mut Vec<i8>,
    columns: &mut impl Iterator<Item = &'a str>,
    format: &str,
    individuals: &[String],
    ploidy: usize,
    num_alleles: usize,
    line: u64,
) -> Result<()> {
    let Some(gt_index) = format.split(':').position(|key| key == "GT") else {
        return Err(Error::VcfDataLine {
            line,
            place: VcfPlace::Column("FORMAT".to_string()),
            problem: format!("`{format}` has no GT key, and GT is the genotype"),
        });
    };
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

impl<R: BufRead + Send> VcfReader<R> {
    /// The next variant of the file, skipping the empty lines and, when the
    /// options ask for it, the variants that failed a filter.
    fn next_variant(&mut self, var: &mut Variant) -> Result<bool> {
        let VcfReader {
            source,
            options,
            individuals,
            chroms,
            needs,
            line,
            line_number,
            finished,
        } = self;
        if *finished {
            return Ok(false);
        }
        loop {
            line.clear();
            if source.read_line(line)? == 0 {
                *finished = true;
                return Ok(false);
            }
            *line_number = next_line_number(*line_number);
            let number = *line_number;
            let text = without_the_line_end(line);
            if text.is_empty() {
                continue;
            }
            // The seven columns up to the FILTER are taken as text and read
            // only when the variant is given: a line that is skipped costs
            // no parsing, and the name of its chromosome gets no number.
            let mut columns = text.split('\t');
            let chrom_text = next_column(&mut columns, "CHROM", number)?;
            let pos_text = next_column(&mut columns, "POS", number)?;
            let id_text = next_column(&mut columns, "ID", number)?;
            let reference_text = next_column(&mut columns, "REF", number)?;
            let alternatives_text = next_column(&mut columns, "ALT", number)?;
            let quality_text = next_column(&mut columns, "QUAL", number)?;
            let filter_text = next_column(&mut columns, "FILTER", number)?;
            if options.only_passed && !passed(filter_text) {
                continue;
            }

            var.pos = parse_position(pos_text, number)?;
            var.chrom = chroms.intern(chrom_text);
            var.filled = Needs::CHROM_POS;
            if needs.contains(Needs::ID) {
                if id_text != MISSING_VALUE {
                    var.id.push_str(id_text);
                }
                var.filled |= Needs::ID;
            }
            // The alleles are counted for every variant that is given, to
            // check the allele numbers of its genotypes, also when the
            // texts of the alleles are not kept.
            let num_alleles = allele_texts(reference_text, alternatives_text).count();
            if needs.contains(Needs::ALLELES) {
                fill_alleles(
                    &mut var.alleles,
                    reference_text,
                    alternatives_text,
                    num_alleles,
                );
                var.filled |= Needs::ALLELES;
            }
            if needs.contains(Needs::QUAL) {
                var.qual = parse_quality(quality_text, number)?;
                var.filled |= Needs::QUAL;
            }
            if needs.contains(Needs::GTS) {
                // INFO is not read, and its column has to be there.
                next_column(&mut columns, "INFO", number)?;
                let format_text = next_column(&mut columns, "FORMAT", number)?;
                fill_genotypes(
                    &mut var.gts,
                    &mut columns,
                    format_text,
                    individuals,
                    options.ploidy,
                    num_alleles,
                    number,
                )?;
                var.filled |= Needs::GTS;
            }
            return Ok(true);
        }
    }
}

impl<R: BufRead + Send> VariantReader for VcfReader<R> {
    fn read_variant(&mut self, var: &mut Variant) -> Result<bool> {
        // Everything but the alleles is emptied here; the alleles are
        // written over in `fill_alleles`, which keeps their strings, and
        // emptied when the reader does not fill them.
        var.chrom = 0;
        var.pos = 0;
        var.gts.clear();
        var.id.clear();
        var.qual = None;
        var.filled = Needs::empty();
        if !self.needs.contains(Needs::ALLELES) {
            var.alleles.clear();
        }
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
        self.needs = needs;
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::{Path, PathBuf};

    use super::{VcfOptions, VcfPlace, VcfReader};
    use crate::error::{Error, Result};
    use crate::variant::{MISSING_ALLELE, Needs, Variant, VariantReader};

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

    /// Every variant a reader gives, until it has no more or it fails.
    fn rows_of(reader: &mut impl VariantReader) -> Result<Vec<Row>> {
        let mut var = Variant::new();
        let mut rows = Vec::new();
        while reader.read_variant(&mut var)? {
            let chrom = reader.chroms().name(var.chrom).unwrap_or("no name");
            rows.push(Row {
                chrom: chrom.to_string(),
                pos: var.pos,
                id: var.id.clone(),
                alleles: var.alleles.clone(),
                qual: var.qual,
                gts: var.gts.clone(),
            });
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
        assert_eq!(
            (line, place),
            (FIRST_DATA_LINE, VcfPlace::Column("FORMAT".to_string()))
        );
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
        assert_eq!(
            (line, place),
            (FIRST_DATA_LINE, VcfPlace::Column("POS".to_string()))
        );
        assert!(problem.contains('x'), "{problem}");
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
