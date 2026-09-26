//! The text of the lines of a VCF that a block holds for the VCF writer of
//! `docs/specs/io_vcf.md`: for each variant, its nine first columns, CHROM
//! to FORMAT, and the column of each individual, as the line has them.
//!
//! Only the VCF reader fills it, and only when it is asked for
//! [`Needs::VCF_TEXT`](crate::variant::Needs::VCF_TEXT). The filters compact
//! it with the other columns of the block and [`Reblock`](super::Reblock)
//! joins and cuts it, so the text of a line stays with the genotypes of its
//! variant: `docs/specs/block.md`.

use std::collections::TryReserveError;

use crate::error::{Error, Result};

/// The text of the line of each variant of a block of a VCF: its nine first
/// columns, CHROM to FORMAT, and the column of each individual, `0/1:12` or
/// `./.`, as the line has them and without its end of line.
///
/// The texts are one buffer for the block, with where each of them ends, and
/// not a string for each individual. The text of a line is the line as it
/// stands after the filters: the filter of individuals takes the columns of
/// the individuals it drops out of it, so it is always the nine first
/// columns and the columns of the individuals of the block, joined by tabs.
#[derive(Debug)]
pub struct VcfText {
    /// How many individuals each line has a column for.
    num_individuals: usize,
    /// The bytes of every line, one after another, each without its end of
    /// line. They are UTF-8: the reader refuses a line that is not.
    /// Compacting the text moves the bytes of a line over the ones of a line
    /// that was dropped, which a `String` does not do without `unsafe`.
    bytes: Vec<u8>,
    /// Where each line ends in `bytes`, one number for each variant.
    line_ends: Vec<usize>,
    /// For each line, one after another, `num_individuals + 1` places
    /// counted from the start of the line: where its nine first columns end,
    /// and where the column of each individual ends. The column of an
    /// individual starts one byte, the tab, after the end before it. They
    /// are counted from the line and not from the block so that compacting
    /// the text moves them as they are, and they are 32 bits so that they
    /// take 4 bytes beside the 4 of a column `0/1` and its tab.
    text_ends: Vec<u32>,
}

impl VcfText {
    /// A text with no line in it, reserved for the ends of `num_vars` lines
    /// of `num_individuals` individuals, and not for their bytes, which
    /// [`VcfText::try_reserve`] asks for a batch of lines at a time.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the ends, which the
    /// reader turns into the error of the crate that names the size of the
    /// block.
    pub(crate) fn with_num_vars(
        num_vars: usize,
        num_individuals: usize,
    ) -> std::result::Result<VcfText, TryReserveError> {
        let mut text = VcfText {
            num_individuals,
            bytes: Vec::new(),
            line_ends: Vec::new(),
            text_ends: Vec::new(),
        };
        text.try_reserve(num_vars, 0)?;
        Ok(text)
    }

    /// How many ends each line has: one for its nine first columns and one
    /// for each individual. The individuals of a block are its genotypes
    /// divided by the ploidy, so they are far below the largest `usize`.
    fn ends_per_line(&self) -> usize {
        self.num_individuals.saturating_add(1)
    }

    /// The memory of `num_lines` more lines holding `num_bytes` bytes of
    /// text, asked of the machine before they are pushed, so that
    /// [`VcfText::push`] does not grow a buffer.
    ///
    /// # Errors
    ///
    /// When the machine does not give it. A number of lines that makes the
    /// ends saturate is far above what any machine gives, and the
    /// reservation of it is the error.
    pub(crate) fn try_reserve(
        &mut self,
        num_lines: usize,
        num_bytes: usize,
    ) -> std::result::Result<(), TryReserveError> {
        self.bytes.try_reserve(num_bytes)?;
        self.line_ends.try_reserve(num_lines)?;
        self.text_ends
            .try_reserve(num_lines.saturating_mul(self.ends_per_line()))?;
        Ok(())
    }

    /// The text of one more line, at the end: `line`, without its end of
    /// line, and `text_ends`, where its nine first columns and the column of
    /// each individual end in it, which the reader found when it parsed the
    /// line.
    ///
    /// # Errors
    ///
    /// When `text_ends` has not one end for the nine first columns and one
    /// for each individual, which is a defect of the reader: the line is
    /// then not added.
    pub(crate) fn push(&mut self, line: &[u8], text_ends: &[u32]) -> Result<()> {
        if text_ends.len() != self.ends_per_line() {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "vcf_text",
                found: text_ends.len(),
                expected: self.ends_per_line(),
            });
        }
        self.bytes.extend_from_slice(line);
        self.line_ends.push(self.bytes.len());
        self.text_ends.extend_from_slice(text_ends);
        Ok(())
    }

    /// Where the line `var` starts and ends in `bytes`, or `None` for a
    /// variant the text does not hold.
    fn bounds_of(&self, var: usize) -> Option<(usize, usize)> {
        let end = self.line_ends.get(var).copied()?;
        let start = match var.checked_sub(1) {
            Some(before) => self.line_ends.get(before).copied()?,
            None => 0,
        };
        Some((start, end))
    }

    /// The bytes of the line `var` and the ends of its texts.
    fn line_of(&self, var: usize) -> Option<(&[u8], &[u32])> {
        let (start, end) = self.bounds_of(var)?;
        let per_line = self.ends_per_line();
        let first = var.checked_mul(per_line)?;
        let ends = self.text_ends.get(first..first.checked_add(per_line)?)?;
        Some((self.bytes.get(start..end)?, ends))
    }

    /// The part of the line `var` from `start` to `end`, as text, and an
    /// empty text when the text holds no such line or no such part.
    fn part_of(&self, var: usize, start: Option<u32>, end: Option<u32>) -> &str {
        let (Some((line, _)), Some(start), Some(end)) = (self.line_of(var), start, end) else {
            return "";
        };
        let (Ok(start), Ok(end)) = (usize::try_from(start), usize::try_from(end)) else {
            return "";
        };
        // The line is UTF-8 and a part begins and ends at a tab or at an end
        // of the line, which are ASCII, so it is UTF-8 too.
        line.get(start..end)
            .and_then(|bytes| std::str::from_utf8(bytes).ok())
            .unwrap_or("")
    }

    /// The ends of the texts of the line `var`, none for a line it does not
    /// hold.
    fn ends_of(&self, var: usize) -> &[u32] {
        self.line_of(var).map_or(&[], |(_, ends)| ends)
    }

    /// How many variants the text holds a line for.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.line_ends.len()
    }

    /// How many individuals each line has a column for.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// CHROM to FORMAT of the variant `var`, counted from 0 in the block,
    /// joined by tabs as the line has them. An empty text for a variant the
    /// block does not hold.
    #[must_use]
    pub fn fixed(&self, var: usize) -> &str {
        self.part_of(var, Some(0), self.ends_of(var).first().copied())
    }

    /// The column of the individual `individual`, counted from 0 in the
    /// block, in the line of the variant `var`, `0/1:12` or `./.`, as the
    /// line has it. An empty text for a variant or an individual the block
    /// does not hold.
    #[must_use]
    pub fn individual(&self, var: usize, individual: usize) -> &str {
        let ends = self.ends_of(var);
        let start = ends
            .get(individual)
            .and_then(|end_before| end_before.checked_add(1));
        let end = individual
            .checked_add(1)
            .and_then(|next| ends.get(next))
            .copied();
        self.part_of(var, start, end)
    }

    /// The columns of every individual of the variant `var`, in the order of
    /// the block, joined by tabs as the line has them, which the writer
    /// copies in one piece after the nine first columns. An empty text for a
    /// variant the block does not hold.
    #[must_use]
    pub fn individuals(&self, var: usize) -> &str {
        let ends = self.ends_of(var);
        let start = ends
            .first()
            .and_then(|end_before| end_before.checked_add(1));
        self.part_of(var, start, ends.last().copied())
    }

    /// The lines whose `keep` is true, in their order, and the rest dropped:
    /// what [`Block::retain_vars`](super::Block::retain_vars) does with the
    /// text. The buffers keep their capacity, and the bytes and the ends of
    /// a line that stays move over the ones of a line that goes.
    ///
    /// `keep` has one value for each line, which the block checked.
    pub(crate) fn retain_vars(&mut self, keep: &[bool]) {
        let per_line = self.ends_per_line();
        let mut write_byte = 0usize;
        let mut write_var = 0usize;
        let mut read_byte = 0usize;
        for (var, keep_it) in keep.iter().enumerate() {
            let Some(line_end) = self.line_ends.get(var).copied() else {
                break;
            };
            if *keep_it {
                // The write position is never past the read one, so no copy
                // writes over bytes or ends that have not been moved yet.
                if read_byte <= line_end && line_end <= self.bytes.len() {
                    self.bytes.copy_within(read_byte..line_end, write_byte);
                }
                write_byte = write_byte.saturating_add(line_end.saturating_sub(read_byte));
                let read_ends = var.saturating_mul(per_line);
                let read_ends_end = read_ends.saturating_add(per_line);
                if read_ends_end <= self.text_ends.len() {
                    self.text_ends
                        .copy_within(read_ends..read_ends_end, write_var.saturating_mul(per_line));
                }
                if let Some(slot) = self.line_ends.get_mut(write_var) {
                    *slot = write_byte;
                }
                write_var = write_var.saturating_add(1);
            }
            read_byte = line_end;
        }
        self.bytes.truncate(write_byte);
        self.line_ends.truncate(write_var);
        self.text_ends.truncate(write_var.saturating_mul(per_line));
    }

    /// The columns of the individuals `keep`, indices into the individuals
    /// of the text, in that order, and the others taken out of every line:
    /// what [`Block::retain_individuals`](super::Block::retain_individuals)
    /// does with the text. Each line becomes its nine first columns and the
    /// kept columns joined by tabs.
    ///
    /// The lines are rewritten one after another, each through a buffer,
    /// since a kept column can come before one that is still to be read. A
    /// line is never longer than it was, so it is written at or before
    /// where it was read and over no line that is still to be read. The two
    /// buffers are allocated once for the block and none for a line.
    ///
    /// `keep` holds each index once, which the block checked.
    ///
    /// # Errors
    ///
    /// When `keep` holds an index at or beyond the individuals of the text,
    /// which the block checked before, and when the ends of the text are
    /// not those of its lines, which the reader and the compactions of this
    /// module build: both are a defect of popnei that no call reaches. The
    /// text is then no longer that of its lines, and the block is lost with
    /// the error.
    pub(crate) fn retain_individuals(&mut self, keep: &[usize]) -> Result<()> {
        let per_line = self.ends_per_line();
        let kept_per_line = keep.len().saturating_add(1);
        // The sizes that the error of a text that is not of its lines names,
        // which do not change until every line has been rewritten.
        let (found, expected) = (
            self.text_ends.len(),
            self.line_ends.len().saturating_mul(per_line),
        );
        let not_of_its_lines = move || Error::BlockArrayOfAnotherSize {
            array: "vcf_text",
            found,
            expected,
        };
        let mut line = Vec::new();
        let mut ends = Vec::with_capacity(kept_per_line);
        let mut write_byte = 0usize;
        let mut read_byte = 0usize;
        for var in 0..self.line_ends.len() {
            let line_end = self
                .line_ends
                .get(var)
                .copied()
                .ok_or_else(not_of_its_lines)?;
            let first_end = var.checked_mul(per_line).ok_or_else(not_of_its_lines)?;
            let read_ends = first_end
                .checked_add(per_line)
                .and_then(|last_end| self.text_ends.get(first_end..last_end))
                .ok_or_else(not_of_its_lines)?;
            let read_line = self
                .bytes
                .get(read_byte..line_end)
                .ok_or_else(not_of_its_lines)?;
            line.clear();
            ends.clear();
            let fixed_end = read_ends.first().copied().ok_or_else(not_of_its_lines)?;
            let fixed = usize::try_from(fixed_end)
                .ok()
                .and_then(|fixed_end| read_line.get(..fixed_end))
                .ok_or_else(not_of_its_lines)?;
            line.extend_from_slice(fixed);
            ends.push(fixed_end);
            for individual in keep {
                let (Some(start), Some(end)) = (
                    read_ends.get(*individual).copied(),
                    individual
                        .checked_add(1)
                        .and_then(|next| read_ends.get(next))
                        .copied(),
                ) else {
                    return Err(Error::IndividualToKeepNotInTheBlock {
                        individual: *individual,
                        num_individuals: self.num_individuals,
                    });
                };
                // The start is the end of the text before the column, the
                // tab, which goes with the column: the kept columns are
                // joined by tabs as the line had them.
                let column = usize::try_from(start)
                    .ok()
                    .zip(usize::try_from(end).ok())
                    .and_then(|(start, end)| read_line.get(start..end))
                    .ok_or_else(not_of_its_lines)?;
                line.extend_from_slice(column);
                // The line is at most as long as the one it came from, whose
                // ends are 32 bits.
                ends.push(u32::try_from(line.len()).map_err(|_| not_of_its_lines())?);
            }
            let write_end = write_byte.saturating_add(line.len());
            self.bytes
                .get_mut(write_byte..write_end)
                .ok_or_else(not_of_its_lines)?
                .copy_from_slice(&line);
            let write_ends = var.saturating_mul(kept_per_line);
            self.text_ends
                .get_mut(write_ends..write_ends.saturating_add(kept_per_line))
                .ok_or_else(not_of_its_lines)?
                .copy_from_slice(&ends);
            if let Some(slot) = self.line_ends.get_mut(var) {
                *slot = write_end;
            }
            write_byte = write_end;
            read_byte = line_end;
        }
        self.bytes.truncate(write_byte);
        self.text_ends
            .truncate(self.line_ends.len().saturating_mul(kept_per_line));
        self.num_individuals = keep.len();
        Ok(())
    }

    /// The lines of `other` after the ones this text holds, which is what
    /// joining two blocks does with their text. Both are of the same
    /// individuals, which [`Block::check`](super::Block::check) said of
    /// each block.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the three buffers.
    pub(crate) fn try_append(
        &mut self,
        other: &VcfText,
    ) -> std::result::Result<(), TryReserveError> {
        // Two buffers that are both allocated have lengths that add without
        // overflow, so this saturating addition never saturates.
        let byte_base = self.bytes.len();
        self.try_reserve(other.line_ends.len(), other.bytes.len())?;
        self.bytes.extend_from_slice(&other.bytes);
        self.line_ends.extend(
            other
                .line_ends
                .iter()
                .map(|end| end.saturating_add(byte_base)),
        );
        self.text_ends.extend_from_slice(&other.text_ends);
        Ok(())
    }

    /// The lines of `count` variants from `from` on, copied into a text of
    /// their own, which is what cutting a block does with them. This text is
    /// left as it is.
    ///
    /// # Errors
    ///
    /// When the machine does not give the memory of the new text.
    pub(crate) fn try_rows(
        &self,
        from: usize,
        count: usize,
    ) -> std::result::Result<VcfText, TryReserveError> {
        let from = from.min(self.line_ends.len());
        let end = from.saturating_add(count).min(self.line_ends.len());
        let byte_from = self
            .bounds_of(from)
            .map_or(self.bytes.len(), |(start, _)| start);
        let mut rows = VcfText::with_num_vars(end.saturating_sub(from), self.num_individuals)?;
        let per_line = self.ends_per_line();
        let bytes = self
            .bytes
            .get(
                byte_from
                    ..end
                        .checked_sub(1)
                        .and_then(|last| self.line_ends.get(last))
                        .copied()
                        .unwrap_or(byte_from),
            )
            .unwrap_or_default();
        rows.bytes.try_reserve_exact(bytes.len())?;
        rows.bytes.extend_from_slice(bytes);
        let ends = self.line_ends.get(from..end).unwrap_or_default();
        rows.line_ends.extend(
            ends.iter()
                .map(|line_end| line_end.saturating_sub(byte_from)),
        );
        let text_ends = self
            .text_ends
            .get(from.saturating_mul(per_line)..end.saturating_mul(per_line))
            .unwrap_or_default();
        rows.text_ends.extend_from_slice(text_ends);
        Ok(rows)
    }
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        reason = "the places of the lines, the variants and the individuals of the small \
                  files of the tests"
    )]

    use std::fs::File;
    use std::io::{BufReader, Cursor};
    use std::path::{Path, PathBuf};

    use crate::block::{Block, BlockReader, Reblock, SourceHeader};
    use crate::error::{Error, Result};
    use crate::filters::{FilteringStats, RegionSelection};
    use crate::io::vcf::{VcfOptions, VcfPlace, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    /// A reference VCF, at the root of the repository beside the Python
    /// tests that read the same files.
    fn reference_vcf(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/vcf")
            .join(name)
    }

    /// The data lines of a reference VCF as the file has them, without
    /// their end of line: what the text of the blocks is checked against,
    /// read with nothing of popnei.
    fn lines_of_the_file(name: &str) -> Vec<String> {
        std::fs::read_to_string(reference_vcf(name))
            .expect("the reference VCF")
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(str::to_owned)
            .collect()
    }

    /// The reader of a reference VCF in blocks of `num_vars_per_block`,
    /// with every line whatever its FILTER, asked for `needs`.
    fn reader_of(
        name: &str,
        num_vars_per_block: usize,
        needs: Needs,
    ) -> VcfReader<BufReader<File>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(num_vars_per_block),
        };
        let mut reader =
            VcfReader::from_path(&reference_vcf(name), options).expect("the reference VCF");
        reader.set_needs(needs);
        reader
    }

    /// Every block a reader gives, each checked.
    fn blocks_of(reader: &mut impl BlockReader) -> Result<Vec<Block>> {
        let mut blocks = Vec::new();
        while let Some(block) = reader.next_block()? {
            block.check()?;
            blocks.push(block);
        }
        Ok(blocks)
    }

    /// The single block of `write.vcf`, every line of it, with every field
    /// and the text.
    fn the_block_of_write_vcf() -> Block {
        let mut reader = reader_of("write.vcf", 100, Needs::ALL | Needs::VCF_TEXT);
        let mut blocks = blocks_of(&mut reader).expect("write.vcf");
        assert_eq!(blocks.len(), 1);
        blocks.pop().expect("the block")
    }

    /// The text of a block, line by line: CHROM to FORMAT and the column of
    /// each individual.
    fn texts_of(block: &Block) -> Vec<(String, Vec<String>)> {
        let text = block.vcf_text.as_ref().expect("the text of the lines");
        (0..text.num_vars())
            .map(|var| {
                let individuals = (0..text.num_individuals())
                    .map(|individual| text.individual(var, individual).to_owned())
                    .collect();
                (text.fixed(var).to_owned(), individuals)
            })
            .collect()
    }

    /// The six lines of `write.vcf`, from "How it is verified" of the writer
    /// in `docs/specs/io_vcf.md`, CHROM to FORMAT and the columns of `a`,
    /// `b` and `c`.
    const THE_SIX_LINES: [(&str, [&str; 3]); 6] = [
        (
            "chr1\t100\trs1\tA\tT\t29.5\tPASS\tAC=4;AN=6;DP=12\tGT:DP",
            ["0/1:4", "0|1:5", "1/1:3"],
        ),
        (
            "chr1\t250\t.\tAT\tA\t.\tq10\tAC=1;AN=4;DP=5\tGT:DP",
            ["./.:0", "0/1:3", "0/0:2"],
        ),
        (
            "chr1\t1000\trs3\tG\tC,T\t50\tPASS\tAC=2,1;AN=6;DP=20\tGT:DP",
            ["1/2:7", "0/1:6", "0/0:7"],
        ),
        (
            "chr1\t1001\t.\tC\t.\t12\t.\tDP=9\tGT",
            ["0/0", "0/0", "0/0"],
        ),
        (
            "chr2\t1\t.\tT\tG\t40\tPASS\tAC=2;AN=5;DP=8\tGT:DP",
            ["1|1:4", "0/.:2", "0/0:2"],
        ),
        (
            "chr2\t1500\trs6\tA\tG\t33\tPASS\tAC=1;AN=6\tGT",
            ["0/0", "1/0", "0/0"],
        ),
    ];

    /// The lines of `THE_SIX_LINES` numbered `lines`, with the columns of
    /// the individuals `individuals` in that order.
    fn the_six_lines(lines: &[usize], individuals: &[usize]) -> Vec<(String, Vec<String>)> {
        lines
            .iter()
            .map(|line| {
                let (fixed, columns) = THE_SIX_LINES[*line];
                let columns = individuals
                    .iter()
                    .map(|individual| columns[*individual].to_owned())
                    .collect();
                (fixed.to_owned(), columns)
            })
            .collect()
    }

    /// That the text of every variant of the block is the line `expected`
    /// with the columns of the individuals `individuals` of the file, and
    /// that it is the line of the variant the other columns of the block
    /// hold at that place: its chromosome, its position and the genotype of
    /// each individual. A text that drifted from its variant fails the
    /// second part although it is a line of the file.
    fn assert_the_text_is_of_its_variants(
        block: &Block,
        chroms: &ChromTable,
        expected: &[&str],
        individuals: &[usize],
    ) {
        let text = block.vcf_text.as_ref().expect("the text of the lines");
        assert_eq!(text.num_vars(), block.num_vars);
        assert_eq!(text.num_vars(), expected.len());
        assert_eq!(text.num_individuals(), individuals.len());
        assert_eq!(block.num_individuals, individuals.len());
        let chrom = block.chrom.as_ref().expect("the chromosomes");
        let pos = block.pos.as_ref().expect("the positions");
        for (var, line) in expected.iter().enumerate() {
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(text.fixed(var), columns[..9].join("\t"), "variant {var}");
            let kept: Vec<&str> = individuals
                .iter()
                .map(|individual| columns[9 + individual])
                .collect();
            assert_eq!(text.individuals(var), kept.join("\t"), "variant {var}");
            assert_eq!(chroms.name(chrom[var]), Some(columns[0]), "variant {var}");
            assert_eq!(pos[var].to_string(), columns[1], "variant {var}");
            let row = &block.gts[var * 2 * individuals.len()..(var + 1) * 2 * individuals.len()];
            for (individual, column) in kept.iter().enumerate() {
                assert_eq!(text.individual(var, individual), *column, "variant {var}");
                let gt = column.split(':').next().expect("the GT");
                let alleles: Vec<i8> = gt
                    .split(['/', '|'])
                    .map(|allele| match allele {
                        "." => MISSING_ALLELE,
                        number => number.parse().expect("an allele"),
                    })
                    .collect();
                assert_eq!(
                    row[individual * 2..individual * 2 + 2],
                    alleles[..],
                    "variant {var}, individual {individual}"
                );
            }
        }
    }

    #[test]
    fn vcf_text_of_write_vcf_is_the_nine_first_columns_and_each_individual_of_its_six_lines() {
        let block = the_block_of_write_vcf();
        assert!(block.fields().contains(Needs::ALL | Needs::VCF_TEXT));
        assert_eq!(
            texts_of(&block),
            the_six_lines(&[0, 1, 2, 3, 4, 5], &[0, 1, 2])
        );
        let text = block.vcf_text.as_ref().expect("the text");
        assert_eq!(text.individuals(0), "0/1:4\t0|1:5\t1/1:3");
        assert_eq!(text.individuals(3), "0/0\t0/0\t0/0");
        // A variant or an individual the block does not hold has an empty
        // text, which no column of a line is.
        assert_eq!(text.fixed(6), "");
        assert_eq!(text.individual(6, 0), "");
        assert_eq!(text.individual(0, 3), "");
        assert_eq!(text.individuals(6), "");
    }

    #[test]
    fn vcf_text_is_in_no_block_that_did_not_ask_for_it_and_alone_in_one_that_asked_for_it_alone() {
        let mut reader = reader_of("write.vcf", 100, Needs::ALL);
        let blocks = blocks_of(&mut reader).expect("write.vcf");
        assert!(blocks[0].vcf_text.is_none());
        assert!(!blocks[0].fields().contains(Needs::VCF_TEXT));

        let mut reader = reader_of("write.vcf", 100, Needs::VCF_TEXT);
        let blocks = blocks_of(&mut reader).expect("write.vcf");
        assert_eq!(blocks[0].fields(), Needs::VCF_TEXT);
        assert!(blocks[0].gts.is_empty());
        assert_eq!(
            texts_of(&blocks[0]),
            the_six_lines(&[0, 1, 2, 3, 4, 5], &[0, 1, 2])
        );
    }

    #[test]
    fn vcf_text_of_write_vcf_read_with_the_default_has_no_line_whose_filter_failed() {
        let options = VcfOptions::default();
        let mut reader =
            VcfReader::from_path(&reference_vcf("write.vcf"), options).expect("write.vcf");
        reader.set_needs(Needs::ALL | Needs::VCF_TEXT);
        let blocks = blocks_of(&mut reader).expect("write.vcf");
        assert_eq!(
            texts_of(&blocks[0]),
            the_six_lines(&[0, 2, 3, 4, 5], &[0, 1, 2])
        );
    }

    /// A VCF of three individuals whose data lines are `lines`, each ended
    /// by `line_end`.
    fn vcf_of(lines: &[&str], line_end: &str) -> Vec<u8> {
        let mut vcf = format!(
            "##fileformat=VCFv4.3{line_end}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc{line_end}"
        );
        for line in lines {
            vcf.push_str(line);
            vcf.push_str(line_end);
        }
        vcf.into_bytes()
    }

    /// The blocks of a VCF held in memory, read with `needs`.
    fn read_with(vcf: Vec<u8>, needs: Needs) -> Result<Vec<Block>> {
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(2),
        };
        let mut reader = VcfReader::new(Cursor::new(vcf), options)?;
        reader.set_needs(needs);
        blocks_of(&mut reader)
    }

    #[test]
    fn vcf_text_of_a_line_that_ends_in_a_carriage_return_has_no_carriage_return() {
        let lines = [
            "chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT:DP\t0/1:3\t./.\t1/1:12",
            "chr1\t9\t.\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1|1",
            "chr2\t1\t.\tC\tG\t.\tPASS\t.\tGT\t./.\t0/1\t0/0",
        ];
        let blocks =
            read_with(vcf_of(&lines, "\r\n"), Needs::GTS | Needs::VCF_TEXT).expect("the VCF");
        assert_eq!(blocks.len(), 2);
        let text = blocks[0].vcf_text.as_ref().expect("the text");
        assert_eq!(text.fixed(0), "chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT:DP");
        assert_eq!(text.individual(0, 2), "1/1:12");
        assert_eq!(text.individual(1, 2), "1|1");
        let text = blocks[1].vcf_text.as_ref().expect("the text");
        assert_eq!(text.individuals(0), "./.\t0/1\t0/0");
    }

    /// The problem of the wrong data line that reading `vcf` with `needs`
    /// gives, with its line and its place.
    fn wrong_line(vcf: Vec<u8>, needs: Needs) -> (u64, VcfPlace, String) {
        match read_with(vcf, needs) {
            Err(Error::VcfDataLine {
                line,
                place,
                problem,
            }) => (line, place, problem),
            other => panic!("not the error of a wrong data line: {other:?}"),
        }
    }

    #[test]
    fn vcf_text_refuses_a_line_that_has_not_one_column_for_each_individual() {
        let good = "chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t./.\t1/1";
        let short = "chr1\t9\t.\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1";
        let long = "chr1\t9\t.\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0\t0/0";
        // Without the genotypes and without the text the columns of the
        // individuals are not counted, and the two lines are read.
        for line in [short, long] {
            let blocks = read_with(vcf_of(&[good, line], "\n"), Needs::CHROM_POS).expect("the VCF");
            assert_eq!(blocks[0].num_vars, 2);
        }
        let (line, place, problem) = wrong_line(vcf_of(&[good, short], "\n"), Needs::VCF_TEXT);
        assert_eq!(line, 4);
        assert_eq!(place, VcfPlace::Line);
        assert_eq!(
            problem,
            "it has the columns of 2 individuals and the header has 3"
        );
        let (line, place, problem) = wrong_line(vcf_of(&[good, long], "\n"), Needs::VCF_TEXT);
        assert_eq!(line, 4);
        assert_eq!(place, VcfPlace::Line);
        assert_eq!(
            problem,
            "it has 2 columns more than the 3 individuals of the header"
        );
    }

    #[test]
    fn vcf_text_refuses_a_column_of_an_individual_whose_bytes_are_not_utf8() {
        let mut vcf = vcf_of(
            &["chr1\t5\t.\tA\tT\t.\tPASS\t.\tGT:XX\t0/1:a\t./.\t1/1"],
            "\n",
        );
        // The `a` of the value of the first individual becomes a byte that
        // is not UTF-8.
        let place = vcf
            .windows(3)
            .position(|bytes| bytes == b"1:a")
            .expect("the value");
        vcf[place + 2] = 0xff;
        // The genotypes are read as bytes, so they alone read the line.
        let blocks = read_with(vcf.clone(), Needs::GTS).expect("the VCF");
        assert_eq!(blocks[0].gts, [0, 1, -1, -1, 1, 1]);
        let (line, place, problem) = wrong_line(vcf, Needs::GTS | Needs::VCF_TEXT);
        assert_eq!(line, 3);
        assert_eq!(place, VcfPlace::Line);
        assert_eq!(problem, "its bytes are not valid UTF-8, and a VCF is text");
    }

    #[test]
    fn vcf_text_retain_vars_keeps_the_lines_of_the_variants_that_stay() {
        let mut block = the_block_of_write_vcf();
        block
            .retain_vars(&[true, false, true, false, false, true])
            .expect("the retain");
        block.check().expect("the block");
        assert_eq!(texts_of(&block), the_six_lines(&[0, 2, 5], &[0, 1, 2]));
        assert_eq!(block.pos, Some(vec![100, 1000, 1500]));

        block
            .retain_vars(&[false, false, true])
            .expect("the retain");
        block.check().expect("the block");
        assert_eq!(texts_of(&block), the_six_lines(&[5], &[0, 1, 2]));

        block.retain_vars(&[false]).expect("the retain");
        block.check().expect("the block");
        assert_eq!(block.num_vars, 0);
        let text = block.vcf_text.as_ref().expect("the text");
        assert_eq!(text.num_vars(), 0);
        assert_eq!(text.num_individuals(), 3);
        assert_eq!(text.fixed(0), "");
    }

    #[test]
    fn vcf_text_retain_vars_that_keeps_every_variant_or_the_last_alone_leaves_their_lines() {
        let mut block = the_block_of_write_vcf();
        block.retain_vars(&[true; 6]).expect("the retain");
        assert_eq!(
            texts_of(&block),
            the_six_lines(&[0, 1, 2, 3, 4, 5], &[0, 1, 2])
        );
        block
            .retain_vars(&[false, false, false, false, false, true])
            .expect("the retain");
        assert_eq!(texts_of(&block), the_six_lines(&[5], &[0, 1, 2]));
        let mut block = the_block_of_write_vcf();
        block
            .retain_vars(&[true, false, false, false, false, false])
            .expect("the retain");
        assert_eq!(texts_of(&block), the_six_lines(&[0], &[0, 1, 2]));
    }

    #[test]
    fn vcf_text_retain_individuals_keeps_the_columns_of_c_and_a_in_that_order() {
        let mut block = the_block_of_write_vcf();
        block.retain_individuals(&[2, 0]).expect("the retain");
        block.check().expect("the block");
        assert_eq!(
            texts_of(&block),
            the_six_lines(&[0, 1, 2, 3, 4, 5], &[2, 0])
        );
        let text = block.vcf_text.as_ref().expect("the text");
        assert_eq!(text.num_individuals(), 2);
        assert_eq!(text.individuals(0), "1/1:3\t0/1:4");
        assert_eq!(text.individuals(5), "0/0\t0/0");
        assert_eq!(text.individual(0, 2), "");
        assert_eq!(block.gts[..4], [1, 1, 0, 1]);

        // A filter of variants after it compacts the lines it left.
        block
            .retain_vars(&[false, true, false, false, true, false])
            .expect("the retain");
        block.check().expect("the block");
        assert_eq!(texts_of(&block), the_six_lines(&[1, 4], &[2, 0]));
    }

    #[test]
    fn vcf_text_retain_individuals_keeps_one_individual_or_every_one_in_another_order() {
        let mut block = the_block_of_write_vcf();
        block.retain_individuals(&[1]).expect("the retain");
        block.check().expect("the block");
        assert_eq!(texts_of(&block), the_six_lines(&[0, 1, 2, 3, 4, 5], &[1]));
        assert_eq!(
            block.vcf_text.as_ref().expect("the text").individuals(4),
            "0/.:2"
        );

        let mut block = the_block_of_write_vcf();
        block.retain_individuals(&[2, 1, 0]).expect("the retain");
        block.check().expect("the block");
        assert_eq!(
            texts_of(&block),
            the_six_lines(&[0, 1, 2, 3, 4, 5], &[2, 1, 0])
        );
    }

    #[test]
    fn vcf_text_retain_individuals_of_a_block_of_no_variants_leaves_a_text_of_the_kept_individuals()
    {
        let mut block = the_block_of_write_vcf();
        block.retain_vars(&[false; 6]).expect("the retain");
        block.retain_individuals(&[2, 0]).expect("the retain");
        block.check().expect("the block");
        let text = block.vcf_text.as_ref().expect("the text");
        assert_eq!(text.num_vars(), 0);
        assert_eq!(text.num_individuals(), 2);
    }

    #[test]
    fn vcf_text_of_a_block_is_checked_against_its_variants_and_its_individuals() {
        let mut block = the_block_of_write_vcf();
        // The text of the block of a filter of individuals that forgot it.
        block.num_individuals = 2;
        block.gts.truncate(6 * 2 * 2);
        match block.check() {
            Err(Error::BlockArrayOfAnotherSize {
                array,
                found,
                expected,
            }) => {
                assert_eq!(array, "individuals of vcf_text");
                assert_eq!((found, expected), (3, 2));
            }
            other => panic!("not the error of a text of other individuals: {other:?}"),
        }

        let mut block = the_block_of_write_vcf();
        // The text of the block of a filter of variants that forgot it.
        block.num_vars = 5;
        block.gts.truncate(5 * 3 * 2);
        block.chrom = None;
        block.pos = None;
        block.id = None;
        block.alleles = None;
        block.qual = None;
        match block.check() {
            Err(Error::BlockArrayOfAnotherSize {
                array,
                found,
                expected,
            }) => {
                assert_eq!(array, "vcf_text");
                assert_eq!((found, expected), (6, 5));
            }
            other => panic!("not the error of a text of other variants: {other:?}"),
        }
    }

    /// The reader over a reader that drops the variants whose place in the
    /// file, counted from 0, `keep` says no to, and then keeps the
    /// individuals `individuals` in that order: what a filter of variants
    /// and a filter of individuals do, with rules whose lines the test
    /// knows. A block it leaves with no variant is not given, as a filter
    /// gives none.
    struct Filtered<R: BlockReader> {
        reader: R,
        keep: fn(usize) -> bool,
        individuals: Vec<usize>,
        names: Vec<String>,
        read: usize,
    }

    impl<R: BlockReader> Filtered<R> {
        fn new(reader: R, keep: fn(usize) -> bool, individuals: &[usize]) -> Filtered<R> {
            let names = individuals
                .iter()
                .map(|individual| reader.individuals()[*individual].clone())
                .collect();
            Filtered {
                reader,
                keep,
                individuals: individuals.to_vec(),
                names,
                read: 0,
            }
        }
    }

    impl<R: BlockReader> BlockReader for Filtered<R> {
        fn next_block(&mut self) -> Result<Option<Block>> {
            while let Some(mut block) = self.reader.next_block()? {
                let keep: Vec<bool> = (self.read..self.read + block.num_vars)
                    .map(self.keep)
                    .collect();
                self.read += block.num_vars;
                block.retain_vars(&keep)?;
                block.retain_individuals(&self.individuals)?;
                if block.num_vars > 0 {
                    return Ok(Some(block));
                }
            }
            Ok(None)
        }

        fn individuals(&self) -> &[String] {
            &self.names
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
            self.reader.filtering_stats()
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

    /// The variants of `many.vcf` that the tests of the filter keep: two of
    /// every three, and none of the third block of 7, the variants 14 to 20,
    /// so that the filter leaves a block with none.
    fn two_of_three_and_not_the_third_block(var: usize) -> bool {
        !var.is_multiple_of(3) && !(14..21).contains(&var)
    }

    /// Where each block of `blocks` starts among the lines `expected`, and
    /// that each has the text of its lines, of the individuals `individuals`.
    fn assert_the_blocks_are_the_lines(
        blocks: &[Block],
        chroms: &ChromTable,
        expected: &[&str],
        individuals: &[usize],
    ) {
        let mut start = 0;
        for block in blocks {
            let end = start + block.num_vars;
            assert_the_text_is_of_its_variants(block, chroms, &expected[start..end], individuals);
            start = end;
        }
        assert_eq!(start, expected.len());
    }

    /// Native alone, since wasm has no threads; the tests of the filter
    /// below read the same blocks of 7 on the one thread there is.
    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn vcf_text_of_many_vcf_in_blocks_of_7_is_each_line_of_the_file_on_one_thread_and_on_several() {
        let lines = lines_of_the_file("many.vcf");
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        assert_eq!(lines.len(), 500);
        let every_individual: Vec<usize> = (0..50).collect();
        for num_threads in [1, 4] {
            for lines_per_batch in [3, 4096] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build()
                    .expect("the pool");
                let mut reader = reader_of("many.vcf", 7, Needs::ALL | Needs::VCF_TEXT);
                reader.set_lines_per_batch(lines_per_batch);
                let blocks = pool.install(|| blocks_of(&mut reader)).expect("many.vcf");
                // 71 blocks of 7 and the last one of 3.
                assert_eq!(blocks.len(), 72);
                assert_eq!(blocks[71].num_vars, 3);
                assert_the_blocks_are_the_lines(
                    &blocks,
                    reader.chroms(),
                    &lines,
                    &every_individual,
                );
            }
        }
    }

    #[test]
    fn vcf_text_of_many_vcf_through_a_filter_and_reblock_stays_with_its_variant() {
        let lines = lines_of_the_file("many.vcf");
        let kept: Vec<&str> = lines
            .iter()
            .enumerate()
            .filter(|(var, _)| two_of_three_and_not_the_third_block(*var))
            .map(|(_, line)| line.as_str())
            .collect();
        assert_eq!(kept.len(), 328);
        // Every individual in the order of the file, and four of them in
        // another order.
        let every_individual: Vec<usize> = (0..50).collect();
        let four_individuals = [49, 3, 17, 0];
        for individuals in [&every_individual[..], &four_individuals[..]] {
            // Blocks of 7 that are joined and cut, of 1, of 5, and one block
            // that holds every variant that is kept, which is the last block
            // and the one that is never cut.
            for (num_vars_per_block, num_blocks, last) in
                [(7, 47, 6), (1, 328, 1), (5, 66, 3), (1000, 1, 328)]
            {
                let reader = Filtered::new(
                    reader_of("many.vcf", 7, Needs::ALL | Needs::VCF_TEXT),
                    two_of_three_and_not_the_third_block,
                    individuals,
                );
                let mut reblock =
                    Reblock::new(reader, Some(num_vars_per_block)).expect("the reblock");
                let blocks = blocks_of(&mut reblock).expect("many.vcf");
                assert_eq!(blocks.len(), num_blocks, "blocks of {num_vars_per_block}");
                assert_eq!(blocks.last().map(|block| block.num_vars), Some(last));
                assert_the_blocks_are_the_lines(&blocks, reblock.chroms(), &kept, individuals);
            }
        }
    }

    #[test]
    fn vcf_text_of_many_vcf_through_a_filter_that_keeps_no_variant_is_no_block() {
        let reader = Filtered::new(
            reader_of("many.vcf", 7, Needs::ALL | Needs::VCF_TEXT),
            |_| false,
            &[1, 0],
        );
        let mut reblock = Reblock::new(reader, Some(7)).expect("the reblock");
        assert_eq!(blocks_of(&mut reblock).expect("many.vcf").len(), 0);
    }

    #[test]
    fn vcf_text_of_many_vcf_with_four_individuals_kept_in_blocks_of_7_is_their_columns() {
        let lines = lines_of_the_file("many.vcf");
        let lines: Vec<&str> = lines.iter().map(String::as_str).collect();
        let individuals = [49, 3, 17, 0];
        let mut reader = reader_of("many.vcf", 7, Needs::ALL | Needs::VCF_TEXT);
        let mut blocks = blocks_of(&mut reader).expect("many.vcf");
        for block in &mut blocks {
            block.retain_individuals(&individuals).expect("the retain");
            block.check().expect("the block");
        }
        assert_the_blocks_are_the_lines(&blocks, reader.chroms(), &lines, &individuals);
    }
}
