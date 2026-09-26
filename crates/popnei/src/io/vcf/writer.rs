//! The VCF writer of `docs/specs/io_vcf.md`: the variants of a pass,
//! after its steps, as a VCF.
//!
//! A variant of a VCF is written as its line was, from the text of the line
//! that the reader kept for it, [`VcfText`]; a variant of a vars file is
//! written from the columns of its block, with what that file holds. The
//! header is the source's, from [`BlockReader::header`], with a `#CHROM`
//! line of the individuals of the pass.

use std::io::Write;

use crate::block::{
    AllelesColumn, Block, BlockReader, GENOTYPES_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK,
    MIN_NUM_VARS_PER_BLOCK, VcfText,
};
use crate::error::{Error, Result};
use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

/// What a line takes out of INFO, and a `##INFO` line out of the header,
/// when the filter of individuals took some out: how often each alternative
/// allele was called and how many alleles were, counts over individuals
/// that are no longer in the file.
const COUNTS_OF_THE_INDIVIDUALS: [&str; 2] = ["AC", "AN"];

/// The nine first columns of the `#CHROM` line.
const CHROM_LINE_COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

/// The header of a VCF written from a source that is not a VCF, before its
/// `##contig` lines, and the one line that comes after them: the format of
/// the file and the one key of the columns of the individuals.
const FILEFORMAT_LINE: &str = "##fileformat=VCFv4.3";
const GT_FORMAT_LINE: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";

/// How many rows of a block one job of the formatting writes at least,
/// into a buffer of its own, which is kept from one block to the next.
///
/// 64 rows of 1000 individuals are 64 lines of 4 bytes an individual, about
/// 256 KB of text in each buffer. Nobody has measured it against another
/// number.
const ROWS_PER_FORMAT_JOB: usize = 64;

mod bgzip;

use bgzip::BgzipOut;

/// How the VCF writer writes its file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VcfWriteOptions {
    /// Members of bgzip, which tabix indexes and bcftools asks a region
    /// of, or plain text.
    pub bgzip: bool,
}

/// What share of the genotypes of a block of the default size,
/// [`GENOTYPES_PER_BLOCK`], a block of a pass of [`write_vcf`] over a VCF
/// holds: a fifth, so that the text of its lines, 20 to 40 bytes for each
/// genotype of a VCF with depths and likelihoods, is 20 to 40 MB. Decided
/// in `docs/specs/io_vcf.md` and not measured.
const SHARE_OF_THE_GENOTYPES_OF_A_BLOCK: usize = 5;

/// The size of the blocks of a pass of [`write_vcf`] over a VCF of
/// `num_individuals` individuals: a fifth of [`GENOTYPES_PER_BLOCK`]
/// divided by the individuals, no fewer than [`MIN_NUM_VARS_PER_BLOCK`] and
/// no more than [`MAX_NUM_VARS_PER_BLOCK`] variants. A source of no
/// individual gives the most, as
/// [`default_num_vars_per_block`](crate::block::default_num_vars_per_block)
/// does.
///
/// A binding crate opens a VCF for a pass of the writer with it, since the
/// block then holds the text of every line beside the genotypes.
#[must_use]
pub fn vcf_text_num_vars_per_block(num_individuals: usize) -> usize {
    (GENOTYPES_PER_BLOCK / SHARE_OF_THE_GENOTYPES_OF_A_BLOCK)
        .checked_div(num_individuals)
        .unwrap_or(MAX_NUM_VARS_PER_BLOCK)
        .clamp(MIN_NUM_VARS_PER_BLOCK, MAX_NUM_VARS_PER_BLOCK)
}

/// Every variant of `reader` into a VCF on `sink`, bgzipped or plain as
/// `options` say, and the sink back with how many variants were written.
///
/// It asks `reader` for every field and for the text of the lines,
/// [`Needs::VCF_TEXT`]. The header is that of `reader.header()`: the lines
/// before `#CHROM` of a VCF, or, from any other source, a `##fileformat`
/// line, one `##contig` line for each chromosome of known length and the
/// `##FORMAT` line of GT; then a `#CHROM` line of `reader.individuals()`.
/// The lines of a VCF are written from the text the reader kept of them;
/// those of any other source from the columns of its blocks, which have to
/// hold the chromosome, the position, the alleles and the genotypes, and an
/// id or a quality the source has no column for is `.`. When
/// the pass has fewer individuals than its source, AC and AN are taken out
/// of every line and their `##INFO` lines out of the header, since they are
/// counts over individuals that are no longer in the file. A source with no
/// variants gives the header alone. Nothing needs a [`Reblock`] before it.
///
/// A bgzipped file is members of 65280 bytes of text, the last one aside,
/// compressed on the threads of rayon and written in order, and the empty
/// member of 28 bytes at the end, which is what tabix indexes and what the
/// reader of popnei refuses a bgzipped file without.
///
/// It borrows the reader, so that the caller reads the counts of the
/// filters of the pass from the chain when it returns.
///
/// [`Reblock`]: crate::block::Reblock
///
/// # Errors
///
/// When the reader fails; when a block is not of its own size or of the
/// individuals and the ploidy of the pass; when a source that is not a VCF
/// has no column of the chromosome, the position, the alleles or the
/// genotypes, [`Error::VcfWriterColumnMissing`]; when a block of a VCF holds
/// no text, or a block a chromosome number its reader has no name for,
/// which are defects of the reader; and when the sink fails, [`Error::VarsFileNotWritten`]. The bytes
/// written before the error are on the sink, and it is the caller that
/// removes the file.
pub fn write_vcf<R: BlockReader + ?Sized, W: Write + Send>(
    reader: &mut R,
    sink: W,
    options: VcfWriteOptions,
) -> Result<(W, u64)> {
    reader.set_needs(Needs::ALL | Needs::VCF_TEXT);
    let num_individuals = reader.individuals().len();
    let ploidy = reader.ploidy();
    let without_counts = num_individuals < reader.header().individuals.len();
    let from = match reader.header().vcf_meta_lines {
        Some(_) => LinesFrom::Text,
        None => LinesFrom::Columns,
    };
    let mut out = match options.bgzip {
        true => VcfOut::Bgzip(BgzipOut::new(sink)),
        false => VcfOut::Plain(sink),
    };
    out.write(&header_of(reader, without_counts))?;
    out.end_of_a_block()?;
    let mut buffers: Vec<Vec<u8>> = Vec::new();
    let mut num_vars: u64 = 0;
    while let Some(block) = reader.next_block()? {
        block.check()?;
        if block.num_individuals != num_individuals || block.ploidy != ploidy {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        let how = LinesOf::block(&block, reader.chroms(), without_counts, from)?;
        format_rows(&how, block.num_vars, &mut buffers)?;
        for buffer in &buffers {
            out.write(buffer)?;
        }
        out.end_of_a_block()?;
        // A variant is a line of the file, so a pass of the
        // 18446744073709551615 variants this count holds is more lines than
        // any file system takes: the sum cannot reach its end. A `usize` is
        // 64 bits natively and 32 in wasm, and both fit in a `u64`.
        num_vars = num_vars.saturating_add(u64::try_from(block.num_vars).unwrap_or(u64::MAX));
    }
    Ok((out.finish()?, num_vars))
}

/// Where the bytes of the file go, in the order of the file: to the sink as
/// they come, or into the members of bgzip.
enum VcfOut<W: Write> {
    Plain(W),
    Bgzip(BgzipOut<W>),
}

impl<W: Write> VcfOut<W> {
    /// `bytes`, after the ones written before.
    ///
    /// # Errors
    ///
    /// When the sink refuses them.
    fn write(&mut self, bytes: &[u8]) -> Result<()> {
        match self {
            VcfOut::Plain(sink) => sink.write_all(bytes).map_err(not_written),
            VcfOut::Bgzip(out) => {
                out.write(bytes);
                Ok(())
            }
        }
    }

    /// The text of the block that was written, compressed into the members
    /// it fills, which is where the threads compress the members of a
    /// bgzipped file; the text that fills no member waits for the next
    /// block.
    ///
    /// # Errors
    ///
    /// When the sink refuses a member.
    fn end_of_a_block(&mut self) -> Result<()> {
        match self {
            VcfOut::Plain(_) => Ok(()),
            VcfOut::Bgzip(out) => out.write_the_full_members(),
        }
    }

    /// The sink, with every byte it was given flushed to it, and in a
    /// bgzipped file the member of the text that was left and the empty
    /// member of the end.
    ///
    /// # Errors
    ///
    /// When the sink refuses them.
    fn finish(self) -> Result<W> {
        match self {
            VcfOut::Plain(mut sink) => {
                sink.flush().map_err(not_written)?;
                Ok(sink)
            }
            VcfOut::Bgzip(out) => out.finish(),
        }
    }
}

/// The error of a file that could not be written, with what the system
/// said and the error it gave, which a binding crate builds its exception
/// with.
fn not_written(failure: std::io::Error) -> Error {
    Error::VarsFileNotWritten {
        problem: failure.to_string(),
        source: Some(failure),
    }
}

/// The header of the file: the lines before `#CHROM`, less the `##INFO`
/// lines of AC and AN when `without_counts`, and the `#CHROM` line of the
/// individuals of the pass, each ended by `\n`.
fn header_of<R: BlockReader + ?Sized>(reader: &R, without_counts: bool) -> Vec<u8> {
    let header = reader.header();
    let mut text = Vec::new();
    match &header.vcf_meta_lines {
        Some(meta_lines) => {
            for line in meta_lines {
                if without_counts && is_the_info_line_of_a_count(line) {
                    continue;
                }
                text.extend_from_slice(line.as_bytes());
                text.push(b'\n');
            }
        }
        None => {
            text.extend_from_slice(FILEFORMAT_LINE.as_bytes());
            text.push(b'\n');
            for (chrom, length) in &header.chrom_lengths {
                text.extend_from_slice(
                    format!("##contig=<ID={chrom},length={length}>\n").as_bytes(),
                );
            }
            text.extend_from_slice(GT_FORMAT_LINE.as_bytes());
            text.push(b'\n');
        }
    }
    text.extend_from_slice(CHROM_LINE_COLUMNS.as_bytes());
    for individual in reader.individuals() {
        text.push(b'\t');
        text.extend_from_slice(individual.as_bytes());
    }
    text.push(b'\n');
    text
}

/// Whether the meta line is the `##INFO` line of AC or of AN: its fields,
/// read as the reader reads those of a `##contig` line, hold an `ID` of one
/// of them. A line that cannot be read that way is not one.
fn is_the_info_line_of_a_count(line: &str) -> bool {
    let Some(fields) = line
        .strip_prefix("##INFO=")
        .map(str::trim)
        .and_then(|value| value.strip_prefix('<'))
        .and_then(|value| value.strip_suffix('>'))
    else {
        return false;
    };
    let Some(fields) = super::fields_outside_quotes(fields) else {
        return false;
    };
    fields.iter().any(|field| {
        field.split_once('=').is_some_and(|(key, value)| {
            key.trim() == "ID" && COUNTS_OF_THE_INDIVIDUALS.contains(&value.trim())
        })
    })
}

/// Where the lines of a pass come from, which the header of its source
/// says: a VCF, whose lines are written from their text, or a source of
/// columns alone, a vars file.
#[derive(Clone, Copy, PartialEq, Eq)]
enum LinesFrom {
    Text,
    Columns,
}

/// What the lines of one block are written from.
enum LinesOf<'a> {
    /// The text the reader kept of each line, and whether AC and AN are
    /// taken out of its INFO.
    Text {
        text: &'a VcfText,
        without_counts: bool,
    },
    /// The columns of the block, each one it has, and the names of its
    /// chromosomes.
    Columns(ColumnsOfABlock<'a>),
}

/// The columns of a block that a line is written from when it has no text:
/// the four every line needs, and the id and the quality, which a source
/// can lack and which are `.` then.
struct ColumnsOfABlock<'a> {
    chroms: &'a ChromTable,
    chrom: &'a [u32],
    pos: &'a [u64],
    id: Option<&'a [String]>,
    alleles: &'a AllelesColumn,
    qual: Option<&'a [f32]>,
    gts: &'a [i8],
    ploidy: usize,
    alleles_per_var: usize,
}

impl<'a> LinesOf<'a> {
    /// The text of the block when it has one, and its columns when it does
    /// not.
    ///
    /// # Errors
    ///
    /// When a block of a VCF, `from` [`LinesFrom::Text`], has no text,
    /// which is a defect of the reader, [`Error::VcfWriterFieldsMissing`];
    /// and when a block of another source lacks the chromosome, the
    /// position, the alleles or the genotypes,
    /// [`Error::VcfWriterColumnMissing`], a source without that column.
    fn block(
        block: &'a Block,
        chroms: &'a ChromTable,
        without_counts: bool,
        from: LinesFrom,
    ) -> Result<Self> {
        if let Some(text) = block.vcf_text.as_ref() {
            return Ok(LinesOf::Text {
                text,
                without_counts,
            });
        }
        if from == LinesFrom::Text {
            return Err(Error::VcfWriterFieldsMissing {
                fields: Needs::VCF_TEXT,
            });
        }
        let missing = |column: &'static str| Error::VcfWriterColumnMissing { column };
        // A block of variants whose `gts` is empty holds no genotypes;
        // `check` said that one that holds them holds all of them.
        if block.gts.is_empty() && block.num_vars > 0 {
            return Err(missing("gts"));
        }
        Ok(LinesOf::Columns(ColumnsOfABlock {
            chroms,
            chrom: block.chrom.as_deref().ok_or_else(|| missing("chrom"))?,
            pos: block.pos.as_deref().ok_or_else(|| missing("pos"))?,
            id: block.id.as_deref(),
            alleles: block.alleles.as_ref().ok_or_else(|| missing("alleles"))?,
            qual: block.qual.as_deref(),
            gts: &block.gts,
            ploidy: block.ploidy.max(1),
            alleles_per_var: block.num_individuals.saturating_mul(block.ploidy),
        }))
    }

    /// The line of the variant `var` of the block, ended by `\n`, after
    /// what `out` holds.
    ///
    /// # Errors
    ///
    /// What [`ColumnsOfABlock::write_line`] refuses.
    fn write_line(&self, var: usize, out: &mut Vec<u8>) -> Result<()> {
        match self {
            LinesOf::Text {
                text,
                without_counts,
            } => {
                let fixed = text.fixed(var);
                match without_counts {
                    true => write_fixed_without_counts(fixed, out),
                    false => out.extend_from_slice(fixed.as_bytes()),
                }
                out.push(b'\t');
                out.extend_from_slice(text.individuals(var).as_bytes());
            }
            LinesOf::Columns(columns) => columns.write_line(var, out)?,
        }
        out.push(b'\n');
        Ok(())
    }
}

/// CHROM to FORMAT of a line, with the values of INFO whose key is AC or AN
/// and the empty ones taken out and `.` for an INFO left with nothing, as `bcftools annotate -x
/// INFO/AC,INFO/AN` writes it.
fn write_fixed_without_counts(fixed: &str, out: &mut Vec<u8>) {
    for (column, text) in fixed.split('\t').enumerate() {
        if column > 0 {
            out.push(b'\t');
        }
        // INFO is the eighth column.
        if column != 7 {
            out.extend_from_slice(text.as_bytes());
            continue;
        }
        let mut written = false;
        for value in text.split(';') {
            let key = value.split_once('=').map_or(value, |(key, _)| key);
            // An empty value, of two `;` side by side or of a `;` at the
            // end, goes as well, as bcftools 1.24 drops it.
            if value.is_empty() || COUNTS_OF_THE_INDIVIDUALS.contains(&key) {
                continue;
            }
            if written {
                out.push(b';');
            }
            out.extend_from_slice(value.as_bytes());
            written = true;
        }
        if !written {
            out.push(b'.');
        }
    }
}

impl ColumnsOfABlock<'_> {
    /// The line of the variant `var`, without its end, after what `out`
    /// holds: the table of "What it gives" of the writer in
    /// `docs/specs/io_vcf.md`.
    ///
    /// # Errors
    ///
    /// [`Error::VcfWriterChromNameMissing`] for a chromosome number the
    /// table has no name for, and [`Error::BlockArrayOfAnotherSize`] for a
    /// column that holds no entry of `var`, which [`Block::check`] made
    /// impossible for a variant of the block.
    fn write_line(&self, var: usize, out: &mut Vec<u8>) -> Result<()> {
        let short = |array: &'static str, found: usize| Error::BlockArrayOfAnotherSize {
            array,
            found,
            expected: var.saturating_add(1),
        };
        let number = *self
            .chrom
            .get(var)
            .ok_or_else(|| short("chrom", self.chrom.len()))?;
        let chrom = self
            .chroms
            .name(number)
            .ok_or(Error::VcfWriterChromNameMissing { number })?;
        out.extend_from_slice(chrom.as_bytes());
        out.push(b'\t');
        let pos = *self
            .pos
            .get(var)
            .ok_or_else(|| short("pos", self.pos.len()))?;
        // A `Vec` takes every byte written into it, so these `write!` cannot
        // fail, and what they format allocates nothing.
        write!(out, "{pos}").map_err(not_written)?;
        out.push(b'\t');
        match self.id.and_then(|id| id.get(var)) {
            Some(id) if !id.is_empty() => out.extend_from_slice(id.as_bytes()),
            _ => out.push(b'.'),
        }
        out.push(b'\t');
        let num_alleles = self.alleles.num_alleles(var);
        if num_alleles == 0 {
            return Err(short("alleles", self.alleles.num_vars()));
        }
        out.extend_from_slice(self.alleles.allele(var, 0).as_bytes());
        out.push(b'\t');
        if num_alleles < 2 {
            out.push(b'.');
        }
        for allele in 1..num_alleles {
            if allele > 1 {
                out.push(b',');
            }
            out.extend_from_slice(self.alleles.allele(var, allele).as_bytes());
        }
        out.push(b'\t');
        match self.qual.and_then(|qual| qual.get(var)) {
            // `Display` of an `f32` is the shortest decimal text that reads
            // back as the same `f32`, and never in the notation with an
            // exponent: `29.5` and `50`.
            Some(qual) if !qual.is_nan() => write!(out, "{qual}").map_err(not_written)?,
            _ => out.push(b'.'),
        }
        out.extend_from_slice(b"\t.\t.\tGT");
        let start = var.saturating_mul(self.alleles_per_var);
        let row = self
            .gts
            .get(start..start.saturating_add(self.alleles_per_var))
            .ok_or_else(|| short("gts", self.gts.len()))?;
        for genotype in row.chunks(self.ploidy) {
            out.push(b'\t');
            for (index, allele) in genotype.iter().enumerate() {
                if index > 0 {
                    out.push(b'/');
                }
                match *allele {
                    MISSING_ALLELE => out.push(b'.'),
                    allele => write!(out, "{allele}").map_err(not_written)?,
                }
            }
        }
        Ok(())
    }
}

/// The lines of the `num_vars` variants of a block, into `buffers`, one
/// buffer for each run of [`ROWS_PER_FORMAT_JOB`] rows, in the order of the
/// rows: the buffers are written one after another. The buffers of the
/// block before are written over, and the ones this block does not need are
/// taken off the end.
///
/// Natively the runs are formatted on the threads of rayon, each into its
/// own buffer, so what is written does not depend on the threads.
///
/// # Errors
///
/// What [`LinesOf::write_line`] refuses, the error of the first run that
/// has one.
fn format_rows(how: &LinesOf<'_>, num_vars: usize, buffers: &mut Vec<Vec<u8>>) -> Result<()> {
    let num_runs = num_vars.div_ceil(ROWS_PER_FORMAT_JOB);
    buffers.truncate(num_runs);
    buffers.resize_with(num_runs, Vec::new);
    let format_run = |(run, buffer): (usize, &mut Vec<u8>)| -> Result<()> {
        buffer.clear();
        let start = run.saturating_mul(ROWS_PER_FORMAT_JOB);
        let end = start.saturating_add(ROWS_PER_FORMAT_JOB).min(num_vars);
        for var in start..end {
            how.write_line(var, buffer)?;
        }
        Ok(())
    };
    #[cfg(not(target_family = "wasm"))]
    {
        use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
        buffers.par_iter_mut().enumerate().try_for_each(format_run)
    }
    #[cfg(target_family = "wasm")]
    {
        buffers.iter_mut().enumerate().try_for_each(format_run)
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        reason = "the sizes and the places of the members of the small files of the tests"
    )]

    use std::fs::File;
    use std::io::{BufReader, Cursor, Write};
    use std::path::{Path, PathBuf};

    use std::process::Command;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use flate2::read::MultiGzDecoder;

    use super::bgzip::TEXT_OF_A_MEMBER;

    /// The empty member that ends a file bgzip wrote, as the bytes of
    /// `printf '' | bgzip -c` of bgzip 1.24, and not the constant of the
    /// writer, which a test would then compare with itself.
    const THE_EMPTY_MEMBER: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    use super::{VcfWriteOptions, vcf_text_num_vars_per_block, write_vcf};
    use crate::block::{Block, BlockReader, SourceHeader};
    use crate::error::{Error, Result};
    use crate::filters::{
        FilteringStats, PassStep, RegionSelection, VarFilteringCriterion, chain_of,
    };
    use crate::io::vars::{VarsReader, write_vars};
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, Needs};

    /// A reference VCF, or what bcftools 1.24 wrote or printed of one, at
    /// the root of the repository.
    fn reference_vcf(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/vcf")
            .join(name)
    }

    /// The text of one of the reference files.
    fn text_of(name: &str) -> String {
        std::fs::read_to_string(reference_vcf(name)).expect("the reference file")
    }

    /// The reader of a reference VCF with `only_passed`, in blocks of
    /// `num_vars_per_block` or of the default size.
    fn reader_of(
        name: &str,
        only_passed: bool,
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed,
            num_vars_per_block,
        };
        VcfReader::from_path(&reference_vcf(name), options).expect("the reference VCF")
    }

    const PLAIN: VcfWriteOptions = VcfWriteOptions { bgzip: false };
    const BGZIP: VcfWriteOptions = VcfWriteOptions { bgzip: true };

    /// The text `write_vcf` writes of `reader` as plain text, and how many
    /// variants it says it wrote.
    fn plain(reader: &mut dyn BlockReader) -> (String, u64) {
        let (bytes, num_vars) = write_vcf(reader, Vec::new(), PLAIN).expect("the VCF was written");
        (String::from_utf8(bytes).expect("text"), num_vars)
    }

    /// The bytes `write_vcf` writes of `reader` bgzipped, and how many
    /// variants it says it wrote.
    fn bgzipped(reader: &mut dyn BlockReader) -> (Vec<u8>, u64) {
        write_vcf(reader, Vec::new(), BGZIP).expect("the VCF was written")
    }

    /// The text `write_vcf` writes of a reader that `make` builds, plain,
    /// and how many variants it says it wrote, after it wrote the same
    /// variants of a second reader of `make` bgzipped and checked that file:
    /// its members, the text they decompress to, which is the plain text,
    /// and what bgzip and tabix say of it.
    fn written(mut make: impl FnMut() -> Box<dyn BlockReader>) -> (String, u64) {
        let (text, num_vars) = plain(&mut *make());
        let (bytes, num_bgzipped) = bgzipped(&mut *make());
        assert_eq!(num_bgzipped, num_vars);
        assert_the_members_are_of_bgzip(&bytes);
        assert_eq!(decompressed(&bytes), text);
        assert_what_bgzip_and_tabix_say(&bytes, &text);
        (text, num_vars)
    }

    /// The text of the members of a bgzipped file, decompressed by flate2
    /// as a gzip stream of many members, which knows nothing of bgzip.
    fn decompressed(bytes: &[u8]) -> String {
        let mut text = String::new();
        std::io::Read::read_to_string(&mut MultiGzDecoder::new(bytes), &mut text)
            .expect("the members decompress");
        text
    }

    /// The members of a bgzipped file, one after another: the bytes of
    /// each and the length of its text, which its last four bytes state.
    /// Each is checked to have the header bgzip writes, with the size of
    /// the member in its extra field `BC`.
    fn members_of(bytes: &[u8]) -> Vec<(&[u8], u32)> {
        let mut members = Vec::new();
        let mut rest = bytes;
        while !rest.is_empty() {
            assert_eq!(rest[..4], [0x1f, 0x8b, 0x08, 0x04], "the start of a member");
            assert_eq!(
                rest[10..16],
                [0x06, 0x00, b'B', b'C', 0x02, 0x00],
                "the extra field"
            );
            let size = usize::from(u16::from_le_bytes([rest[16], rest[17]])) + 1;
            let (member, after) = rest.split_at(size);
            let length = u32::from_le_bytes(member[size - 4..].try_into().expect("four bytes"));
            members.push((member, length));
            rest = after;
        }
        members
    }

    /// That every member of a bgzipped file but the last two holds 65280
    /// bytes of text, that the one before the last holds the rest, and
    /// that the last is the empty member of 28 bytes that htslib writes.
    fn assert_the_members_are_of_bgzip(bytes: &[u8]) {
        let members = members_of(bytes);
        let (last, texts) = members.split_last().expect("a member");
        assert_eq!(last.0, THE_EMPTY_MEMBER);
        assert_eq!(last.0.len(), 28);
        let (rest, full) = texts.split_last().expect("a member of text");
        assert!(rest.1 > 0 && rest.1 as usize <= TEXT_OF_A_MEMBER);
        for (_, length) in full {
            assert_eq!(*length as usize, TEXT_OF_A_MEMBER);
        }
    }

    /// Which of the files of these tests this is, so that two tests that
    /// run at the same time write two files.
    static FILES_WRITTEN: AtomicUsize = AtomicUsize::new(0);

    /// Whether `program` can be run, which is how a test that runs it says
    /// that it is skipped: cargo has no way to mark a test skipped, so the
    /// test passes and prints the program that was not there.
    fn can_run(program: &str) -> bool {
        let found = Command::new(program).arg("--version").output().is_ok();
        if !found {
            eprintln!("skipped: {program} is not in the PATH");
        }
        found
    }

    /// That `bgzip -t` finds the bgzipped file whole, that `tabix -p vcf`
    /// indexes it, and that `tabix` on `chr1:900-1100` gives the lines of
    /// the plain `text` on chr1 from 900 to 1100, each with its end. The
    /// tests that write `write.vcf` and `many.vcf` check those lines
    /// against the files.
    fn assert_what_bgzip_and_tabix_say(bytes: &[u8], text: &str) {
        if !can_run("bgzip") || !can_run("tabix") {
            return;
        }
        let path = std::env::temp_dir().join(format!(
            "popnei-write-vcf-{}-{}.vcf.gz",
            std::process::id(),
            FILES_WRITTEN.fetch_add(1, Ordering::SeqCst)
        ));
        let index = path.with_extension("gz.tbi");
        std::fs::write(&path, bytes).expect("the file was written");
        let file = path.to_str().expect("a path of text");
        let run = |program: &str, args: &[&str]| {
            let output = Command::new(program)
                .args(args)
                .output()
                .expect("the program ran");
            assert!(output.status.success(), "{program} {args:?}: {output:?}");
            String::from_utf8(output.stdout).expect("text")
        };
        run("bgzip", &["-t", file]);
        run("tabix", &["-f", "-p", "vcf", file]);
        let region = run("tabix", &[file, "chr1:900-1100"]);
        let expected: String = text
            .split_inclusive('\n')
            .filter(|line| {
                let mut columns = line.split('\t');
                let chrom = columns.next();
                let pos = columns.next().and_then(|pos| pos.parse::<u64>().ok());
                chrom == Some("chr1") && pos.is_some_and(|pos| (900..=1100).contains(&pos))
            })
            .collect();
        assert_eq!(region, expected);
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_file(&index);
    }

    /// The text of `write.vcf` of "How it is verified" of the writer in
    /// `docs/specs/io_vcf.md`, which the tests write back.
    fn write_vcf_text() -> String {
        text_of("write.vcf")
    }

    #[test]
    fn write_vcf_of_write_vcf_read_with_every_line_gives_its_bytes_in_blocks_of_any_size() {
        for num_vars_per_block in [None, Some(1), Some(4)] {
            let (text, num_vars) =
                written(|| Box::new(reader_of("write.vcf", false, num_vars_per_block)));
            assert_eq!(text, write_vcf_text(), "blocks of {num_vars_per_block:?}");
            assert_eq!(num_vars, 6);
        }
    }

    #[test]
    fn write_vcf_of_many_vcf_in_blocks_of_7_and_of_the_default_gives_its_bytes_on_any_threads() {
        let many = text_of("many.vcf");
        let mut one_thread_bytes: Option<Vec<u8>> = None;
        for num_vars_per_block in [Some(7), None] {
            for num_threads in [1, 4] {
                #[cfg(not(target_family = "wasm"))]
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build()
                    .expect("the pool");
                let make = || -> Box<dyn BlockReader> {
                    Box::new(reader_of("many.vcf", false, num_vars_per_block))
                };
                #[cfg(not(target_family = "wasm"))]
                let ((text, num_vars), bytes) =
                    pool.install(|| (written(make), bgzipped(&mut *make()).0));
                #[cfg(target_family = "wasm")]
                let ((text, num_vars), bytes) = {
                    let _ = num_threads;
                    (written(make), bgzipped(&mut *make()).0)
                };
                assert_eq!(text, many, "blocks of {num_vars_per_block:?}");
                assert_eq!(num_vars, 500);
                // The bytes of the members do not depend on the threads
                // that compressed them.
                let first = one_thread_bytes.get_or_insert_with(|| bytes.clone());
                assert_eq!(*first, bytes);
            }
        }
    }

    #[test]
    fn write_vcf_of_write_vcf_read_with_the_default_leaves_out_the_line_whose_filter_failed() {
        let (text, num_vars) = written(|| Box::new(reader_of("write.vcf", true, None)));
        let expected: String = write_vcf_text()
            .split_inclusive('\n')
            .filter(|line| !line.starts_with("chr1\t250\t"))
            .collect();
        assert_eq!(expected.lines().count(), 15);
        assert_eq!(text, expected);
        assert_eq!(num_vars, 5);
    }

    /// The text `write_vcf` writes of `write.vcf` read with every line and
    /// with the filter of individuals that keeps `names`.
    fn written_with_the_individuals(names: &[&str]) -> String {
        let names: Vec<String> = names.iter().map(|name| (*name).to_owned()).collect();
        written(|| {
            let reader = reader_of("write.vcf", false, None);
            chain_of(
                Box::new(reader),
                &[PassStep::KeepIndividuals(names.clone())],
            )
            .expect("the chain")
        })
        .0
    }

    #[test]
    fn write_vcf_with_the_filter_of_c_and_a_is_what_bcftools_writes_less_its_pass_line() {
        let text = written_with_the_individuals(&["c", "a"]);
        let bcftools = text_of("write.c_a.bcftools.vcf");
        let pass_line = "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
        assert_eq!(bcftools.split_inclusive('\n').nth(1), Some(pass_line));
        assert_eq!(text, bcftools.replacen(pass_line, "", 1));
        assert!(text.contains("\nchr1\t100\trs1\tA\tT\t29.5\tPASS\tDP=12\tGT:DP\t1/1:3\t0/1:4\n"));
        assert!(text.ends_with("\nchr2\t1500\trs6\tA\tG\t33\tPASS\t.\tGT\t0/0\t0/0\n"));
    }

    #[test]
    fn write_vcf_with_the_filter_of_c_b_and_a_keeps_ac_and_an() {
        let text = written_with_the_individuals(&["c", "b", "a"]);
        // write.vcf with the three columns of the individuals of every line
        // after the header in the other order.
        let expected: String = write_vcf_text()
            .lines()
            .map(|line| {
                if line.starts_with("##") {
                    return format!("{line}\n");
                }
                let mut columns: Vec<&str> = line.split('\t').collect();
                columns[9..].reverse();
                format!("{}\n", columns.join("\t"))
            })
            .collect();
        assert!(expected.contains("\tAC=4;AN=6;DP=12\tGT:DP\t1/1:3\t0|1:5\t0/1:4\n"));
        assert_eq!(text, expected);
    }

    #[test]
    fn write_vcf_from_the_vars_file_of_write_vcf_gives_the_five_lines_of_the_spec() {
        let reader = reader_of("write.vcf", true, None);
        let (vars, _) = write_vars(reader, Vec::new(), None).expect("the vars file");
        let (text, num_vars) = written(|| {
            Box::new(VarsReader::new(Cursor::new(vars.clone())).expect("the vars file"))
        });
        let expected = "##fileformat=VCFv4.3\n\
            ##contig=<ID=chr1,length=2000>\n\
            ##contig=<ID=chr2,length=1500>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n\
            chr1\t100\trs1\tA\tT\t29.5\t.\t.\tGT\t0/1\t0/1\t1/1\n\
            chr1\t1000\trs3\tG\tC,T\t50\t.\t.\tGT\t1/2\t0/1\t0/0\n\
            chr1\t1001\t.\tC\t.\t12\t.\t.\tGT\t0/0\t0/0\t0/0\n\
            chr2\t1\t.\tT\tG\t40\t.\t.\tGT\t1/1\t0/.\t0/0\n\
            chr2\t1500\trs6\tA\tG\t33\t.\t.\tGT\t0/0\t1/0\t0/0\n";
        assert_eq!(text, expected);
        assert_eq!(num_vars, 5);
        // The rows `bcftools query` prints of write.vcf read with the
        // default, but for the two separators of the phase, which a vars
        // file does not keep.
        let rows: Vec<String> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let columns: Vec<&str> = line.split('\t').collect();
                let mut row = columns[..6].to_vec();
                row.extend(&columns[9..]);
                row.join("\t")
            })
            .collect();
        let bcftools: Vec<String> = text_of("write.passed.bcftools.tsv")
            .lines()
            .map(|row| row.replace('|', "/"))
            .collect();
        assert_eq!(rows, bcftools);
    }

    #[test]
    fn write_vcf_of_many_vcf_with_the_missing_data_filter_at_0_04_is_the_215_lines_of_bcftools() {
        let (text, num_vars) = written(|| {
            let reader = reader_of("many.vcf", false, Some(7));
            chain_of(
                Box::new(reader),
                &[PassStep::VarFilter(VarFilteringCriterion::MaxMissingRate(
                    0.04,
                ))],
            )
            .expect("the chain")
        });
        let header: String = text_of("many.vcf")
            .split_inclusive('\n')
            .take_while(|line| line.starts_with('#'))
            .collect();
        let bcftools = text_of("many.missing_0.04.bcftools.vcf");
        assert_eq!(bcftools.lines().count(), 215);
        assert_eq!(text, header + &bcftools);
        assert_eq!(num_vars, 215);
    }

    #[test]
    fn write_vcf_of_a_source_with_no_variant_gives_the_header_alone() {
        let header = "##fileformat=VCFv4.3\n\
            ##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n";
        let written = written(|| {
            Box::new(
                VcfReader::new(
                    Cursor::new(header.as_bytes().to_vec()),
                    VcfOptions::default(),
                )
                .expect("the VCF"),
            )
        });
        assert_eq!(written, (header.to_owned(), 0));
    }

    #[test]
    fn write_vcf_takes_out_the_info_lines_of_ac_and_an_by_their_id_and_no_other() {
        let lines = [
            "##fileformat=VCFv4.3",
            "##INFO=< ID=AN ,Number=1,Type=Integer,Description=\"Allele number\">",
            "##INFO=<ID=ACX,Number=1,Type=Integer,Description=\"ID=AC, not AC\">",
            "##INFO=<Number=A,ID=AC,Type=Integer,Description=\"Allele count\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
            "chr1\t5\t.\tA\tT\t.\tPASS\tAC;ACX=2;AN=4;AF=0.5\tGT\t0/1\t1/1",
            "chr1\t6\t.\tA\tT\t.\tPASS\tAN=4\tGT\t0/0\t1/1",
            "chr1\t7\t.\tA\tT\t.\tPASS\tAC=1;AN=4;\tGT\t0/0\t1/1",
            "chr1\t8\t.\tA\tT\t.\tPASS\tAC=1;;DP=3\tGT\t0/0\t1/1",
        ];
        let vcf = lines.join("\n") + "\n";
        let make = || {
            let reader =
                VcfReader::new(Cursor::new(vcf.clone().into_bytes()), VcfOptions::default())
                    .expect("the VCF");
            chain_of(
                Box::new(reader),
                &[PassStep::KeepIndividuals(vec!["b".to_owned()])],
            )
            .expect("the chain")
        };
        let expected = [
            "##fileformat=VCFv4.3",
            "##INFO=<ID=ACX,Number=1,Type=Integer,Description=\"ID=AC, not AC\">",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tb",
            "chr1\t5\t.\tA\tT\t.\tPASS\tACX=2;AF=0.5\tGT\t1/1",
            "chr1\t6\t.\tA\tT\t.\tPASS\t.\tGT\t1/1",
            // bcftools 1.24 writes the INFO of these two as `.` and `DP=3`:
            // the empty values a `;` at the end and two side by side leave
            // go with AC and AN.
            "chr1\t7\t.\tA\tT\t.\tPASS\t.\tGT\t1/1",
            "chr1\t8\t.\tA\tT\t.\tPASS\tDP=3\tGT\t1/1",
        ];
        assert_eq!(written(make).0, expected.join("\n") + "\n");
    }

    /// A reader that gives one block it was built with, over the individuals
    /// and the header of `write.vcf`.
    struct OneBlock {
        block: Option<Block>,
        individuals: Vec<String>,
        chroms: ChromTable,
        header: SourceHeader,
    }

    impl OneBlock {
        /// The first block of `write.vcf` read with `needs`, which the
        /// reader it came from was asked for however `write_vcf` asks.
        fn of_write_vcf(needs: Needs) -> OneBlock {
            let mut reader = reader_of("write.vcf", false, None);
            reader.set_needs(needs);
            let block = reader.next_block().expect("write.vcf");
            OneBlock {
                block,
                individuals: reader.individuals().to_vec(),
                chroms: reader.chroms().clone(),
                header: reader.header().clone(),
            }
        }
    }

    impl BlockReader for OneBlock {
        fn next_block(&mut self) -> Result<Option<Block>> {
            Ok(self.block.take())
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

        fn set_needs(&mut self, _needs: Needs) {}

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }

        fn header(&self) -> &SourceHeader {
            &self.header
        }

        fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
            false
        }

        fn num_skipped(&self) -> u64 {
            0
        }
    }

    /// The reader over `write.vcf`, read with the default, whose blocks lose
    /// the columns of `dropped`: a source of fewer fields, as the vars file
    /// of another writer can be.
    struct WithoutColumns {
        reader: VcfReader<BufReader<File>>,
        dropped: Needs,
    }

    impl BlockReader for WithoutColumns {
        fn next_block(&mut self) -> Result<Option<Block>> {
            let Some(mut block) = self.reader.next_block()? else {
                return Ok(None);
            };
            if self.dropped.contains(Needs::CHROM_POS) {
                block.chrom = None;
                block.pos = None;
            }
            if self.dropped.contains(Needs::ID) {
                block.id = None;
            }
            if self.dropped.contains(Needs::ALLELES) {
                block.alleles = None;
            }
            if self.dropped.contains(Needs::QUAL) {
                block.qual = None;
            }
            Ok(Some(block))
        }

        fn individuals(&self) -> &[String] {
            self.reader.individuals()
        }

        fn ploidy(&self) -> usize {
            self.reader.ploidy()
        }

        fn chroms(&self) -> &ChromTable {
            self.reader.chroms()
        }

        fn set_needs(&mut self, needs: Needs) {
            self.reader.set_needs(needs);
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }

        fn header(&self) -> &SourceHeader {
            self.reader.header()
        }

        fn skip_outside(&mut self, _selection: RegionSelection) -> bool {
            false
        }

        fn num_skipped(&self) -> u64 {
            0
        }
    }

    /// The vars file of `write.vcf` read with the default, written from
    /// blocks without the columns of `dropped`.
    fn vars_file_without(dropped: Needs) -> Vec<u8> {
        let reader = WithoutColumns {
            reader: reader_of("write.vcf", true, None),
            dropped,
        };
        write_vars(reader, Vec::new(), None)
            .expect("the vars file")
            .0
    }

    #[test]
    fn write_vcf_of_a_vars_file_without_the_id_and_the_qual_writes_a_dot_for_both() {
        let vars = vars_file_without(Needs::ID | Needs::QUAL);
        let (text, num_vars) = written(|| {
            Box::new(VarsReader::new(Cursor::new(vars.clone())).expect("the vars file"))
        });
        let lines: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
        assert_eq!(
            lines,
            [
                "chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/1\t1/1",
                "chr1\t1000\t.\tG\tC,T\t.\t.\t.\tGT\t1/2\t0/1\t0/0",
                "chr1\t1001\t.\tC\t.\t.\t.\t.\tGT\t0/0\t0/0\t0/0",
                "chr2\t1\t.\tT\tG\t.\t.\t.\tGT\t1/1\t0/.\t0/0",
                "chr2\t1500\t.\tA\tG\t.\t.\t.\tGT\t0/0\t1/0\t0/0",
            ]
        );
        assert_eq!(num_vars, 5);
    }

    #[test]
    fn write_vcf_of_a_vars_file_without_the_alleles_or_the_chrom_names_the_column() {
        for (dropped, column) in [(Needs::ALLELES, "alleles"), (Needs::CHROM_POS, "chrom")] {
            let vars = vars_file_without(dropped);
            let mut reader = VarsReader::new(Cursor::new(vars)).expect("the vars file");
            match write_vcf(&mut reader, Vec::new(), PLAIN) {
                Err(error @ Error::VcfWriterColumnMissing { .. }) => {
                    assert!(
                        error.to_string().contains(&format!("`{column}`")),
                        "{error}"
                    );
                    assert!(error.names_the_file());
                }
                other => panic!("not the error of a missing {column}: {other:?}"),
            }
        }
    }

    #[test]
    fn write_vcf_writes_a_block_of_a_source_that_is_not_a_vcf_from_its_columns() {
        let text = written(|| {
            let mut reader = OneBlock::of_write_vcf(Needs::ALL);
            // The header of a source that is not a VCF has no lines of one.
            reader.header.vcf_meta_lines = None;
            Box::new(reader)
        })
        .0;
        assert!(text.contains("\nchr1\t100\trs1\tA\tT\t29.5\t.\t.\tGT\t0/1\t0/1\t1/1\n"));
        assert!(text.contains("\nchr1\t250\t.\tAT\tA\t.\t.\t.\tGT\t./.\t0/1\t0/0\n"));
    }

    #[test]
    fn write_vcf_refuses_a_block_of_a_vcf_without_its_text() {
        let mut reader = OneBlock::of_write_vcf(Needs::ALL);
        match write_vcf(&mut reader, Vec::new(), PLAIN) {
            Err(Error::VcfWriterFieldsMissing { fields }) => {
                assert_eq!(fields, Needs::VCF_TEXT);
            }
            other => panic!("not the error of the missing text: {other:?}"),
        }
    }

    #[test]
    fn write_vcf_refuses_a_chromosome_number_its_reader_has_no_name_for() {
        let mut reader = OneBlock::of_write_vcf(Needs::ALL);
        reader.header.vcf_meta_lines = None;
        reader.chroms = ChromTable::new();
        match write_vcf(&mut reader, Vec::new(), PLAIN) {
            Err(Error::VcfWriterChromNameMissing { number }) => assert_eq!(number, 0),
            other => panic!("not the error of a chromosome with no name: {other:?}"),
        }
    }

    #[test]
    fn write_vcf_refuses_a_block_of_other_individuals_than_its_reader_gives() {
        let mut reader = OneBlock::of_write_vcf(Needs::ALL | Needs::VCF_TEXT);
        reader.individuals.pop();
        match write_vcf(&mut reader, Vec::new(), PLAIN) {
            Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                found_num_individuals,
                ..
            }) => assert_eq!((num_individuals, found_num_individuals), (2, 3)),
            other => panic!("not the error of blocks that do not fit: {other:?}"),
        }
    }

    #[test]
    fn write_vcf_from_a_vars_file_writes_no_contig_line_for_a_chromosome_with_no_length() {
        let vcf = "##fileformat=VCFv4.3\n\
            ##contig=<ID=chr1,length=2000>\n\
            ##contig=<ID=chr2>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
            chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n\
            chr2\t7\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n";
        let reader = VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), VcfOptions::default())
            .expect("the VCF");
        let (vars, _) = write_vars(reader, Vec::new(), None).expect("the vars file");
        let (text, _) = written(|| {
            Box::new(VarsReader::new(Cursor::new(vars.clone())).expect("the vars file"))
        });
        let expected = "##fileformat=VCFv4.3\n\
            ##contig=<ID=chr1,length=2000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\n\
            chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/1\n\
            chr2\t7\t.\tA\tT\t.\t.\t.\tGT\t1/1\n";
        assert_eq!(text, expected);
    }

    /// A sink that takes `room` bytes and refuses the rest.
    #[derive(Debug)]
    struct SinkThatFills {
        room: usize,
    }

    impl Write for SinkThatFills {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if self.room == 0 {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::StorageFull,
                    "the disc is full",
                ));
            }
            let taken = bytes.len().min(self.room);
            self.room = self.room.saturating_sub(taken);
            Ok(taken)
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn write_vcf_into_a_sink_that_fills_is_the_error_of_a_file_not_written() {
        for (room, options) in [
            (0, PLAIN),
            (100, PLAIN),
            (800, PLAIN),
            (0, BGZIP),
            (300, BGZIP),
        ] {
            let mut reader = reader_of("write.vcf", false, None);
            match write_vcf(&mut reader, SinkThatFills { room }, options) {
                Err(Error::VarsFileNotWritten { source, .. }) => {
                    let kind = source.map(|source| source.kind());
                    assert_eq!(kind, Some(std::io::ErrorKind::StorageFull));
                }
                other => panic!("not the error of a file not written: {other:?}"),
            }
        }
    }

    #[test]
    fn write_vcf_bgzipped_of_many_vcf_is_a_member_of_65280_bytes_the_rest_and_the_empty_one() {
        let (bytes, _) = bgzipped(&mut reader_of("many.vcf", false, Some(7)));
        let lengths: Vec<u32> = members_of(&bytes)
            .iter()
            .map(|(_, length)| *length)
            .collect();
        // many.vcf is 117346 bytes.
        assert_eq!(lengths, [65280, 52066, 0]);
        assert_eq!(members_of(&bytes)[2].0, THE_EMPTY_MEMBER);
        assert_eq!(decompressed(&bytes), text_of("many.vcf"));
    }

    #[test]
    fn write_vcf_bgzipped_of_the_six_lines_of_write_vcf_is_one_member_and_the_empty_one() {
        let (bytes, _) = bgzipped(&mut reader_of("write.vcf", false, None));
        let lengths: Vec<u32> = members_of(&bytes)
            .iter()
            .map(|(_, length)| *length)
            .collect();
        assert_eq!(lengths, [845, 0]);
        assert!(bytes.ends_with(&THE_EMPTY_MEMBER));
    }

    /// The data lines `write_vcf` writes from the vars file of the VCF of
    /// two individuals whose data lines are `lines`, read with `ploidy`.
    fn lines_from_the_vars_file(ploidy: usize, lines: &[&str]) -> Vec<String> {
        let vcf = format!(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n{}\n",
            lines.join("\n")
        );
        let options = VcfOptions {
            ploidy,
            only_passed: false,
            num_vars_per_block: None,
        };
        let reader = VcfReader::new(Cursor::new(vcf.into_bytes()), options).expect("the VCF");
        let (vars, _) = write_vars(reader, Vec::new(), None).expect("the vars file");
        let (text, _) = written(|| {
            Box::new(VarsReader::new(Cursor::new(vars.clone())).expect("the vars file"))
        });
        text.lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect()
    }

    #[test]
    fn write_vcf_from_a_vars_file_writes_haploid_and_tetraploid_genotypes() {
        let haploid = lines_from_the_vars_file(
            1,
            &[
                "chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT\t0\t1",
                "chr1\t9\t.\tA\tT,G\t.\tPASS\t.\tGT\t2\t.",
            ],
        );
        assert_eq!(
            haploid,
            [
                "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0\t1",
                "chr1\t9\t.\tA\tT,G\t.\t.\t.\tGT\t2\t.",
            ]
        );
        let tetraploid = lines_from_the_vars_file(
            4,
            &[
                "chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT\t0/0/1/1\t0|1|1|1",
                "chr1\t9\t.\tA\tT\t.\tPASS\t.\tGT\t.\t0/./1/.",
            ],
        );
        assert_eq!(
            tetraploid,
            [
                "chr1\t5\t.\tA\tT\t.\t.\t.\tGT\t0/0/1/1\t0/1/1/1",
                "chr1\t9\t.\tA\tT\t.\t.\t.\tGT\t./././.\t0/./1/.",
            ]
        );
    }

    /// A VCF of two individuals and `num_lines` data lines, whose bytes are
    /// `num_bytes` when that is given: the id of its last line is made as
    /// long as it takes.
    fn vcf_of_lines(num_lines: usize, num_bytes: Option<usize>) -> String {
        let mut vcf = String::from(
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n",
        );
        for line in 0..num_lines {
            let id = format!("var{line}");
            let id = match num_bytes {
                Some(num_bytes) if line + 1 == num_lines => {
                    let rest = format!("chr1\t{}\t\tA\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n", line + 1);
                    "x".repeat(num_bytes - vcf.len() - rest.len())
                }
                _ => id,
            };
            vcf.push_str(&format!(
                "chr1\t{}\t{id}\tA\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n",
                line + 1
            ));
        }
        vcf
    }

    /// The bytes `write_vcf` writes of `vcf`, bgzipped, on `num_threads`.
    #[cfg(not(target_family = "wasm"))]
    fn bgzipped_on(vcf: &str, num_threads: usize) -> Vec<u8> {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .expect("the pool");
        pool.install(|| {
            let mut reader =
                VcfReader::new(Cursor::new(vcf.as_bytes().to_vec()), VcfOptions::default())
                    .expect("the VCF");
            bgzipped(&mut reader).0
        })
    }

    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn write_vcf_bgzipped_of_a_block_of_five_members_writes_them_in_order_on_any_threads() {
        // 7000 lines of 40 to 44 bytes, in one block of the default size,
        // 10000 variants for two individuals: five members of text.
        let vcf = vcf_of_lines(7000, None);
        assert!(
            vcf.len() > 4 * 65280 && vcf.len() < 5 * 65280,
            "{}",
            vcf.len()
        );
        let one_thread = bgzipped_on(&vcf, 1);
        let lengths: Vec<u32> = members_of(&one_thread)
            .iter()
            .map(|(_, length)| *length)
            .collect();
        assert_eq!(lengths.len(), 6);
        assert_eq!(decompressed(&one_thread), vcf);
        assert_eq!(bgzipped_on(&vcf, 4), one_thread);
        assert_the_members_are_of_bgzip(&one_thread);
    }

    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn write_vcf_bgzipped_of_a_text_of_two_members_exactly_ends_with_no_member_of_the_rest() {
        let vcf = vcf_of_lines(3000, Some(2 * 65280));
        assert_eq!(vcf.len(), 2 * 65280);
        let bytes = bgzipped_on(&vcf, 4);
        let lengths: Vec<u32> = members_of(&bytes)
            .iter()
            .map(|(_, length)| *length)
            .collect();
        assert_eq!(lengths, [65280, 65280, 0]);
        assert_eq!(decompressed(&bytes), vcf);
    }

    #[test]
    fn vcf_text_num_vars_per_block_is_a_fifth_of_the_genotypes_of_a_block_from_100_to_10000() {
        assert_eq!(vcf_text_num_vars_per_block(1000), 1000);
        assert_eq!(vcf_text_num_vars_per_block(3000), 333);
        assert_eq!(vcf_text_num_vars_per_block(50), 10000);
        assert_eq!(vcf_text_num_vars_per_block(20000), 100);
        assert_eq!(vcf_text_num_vars_per_block(0), 10000);
    }
}
