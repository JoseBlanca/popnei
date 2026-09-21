//! The error of the crate.
//!
//! Every module of popnei adds its cases to [`Error`], and every operation
//! that can fail returns the [`Result`] of this module. The enum is
//! `#[non_exhaustive]`, so a module that is written later adds a case
//! without breaking the code that matches on it.

use std::path::PathBuf;

use thiserror::Error as ThisError;

use crate::block::BlockSize;
use crate::io::vcf::VcfPlace;
use crate::variant::Needs;

/// Anything that went wrong in popnei.
#[derive(Debug, ThisError)]
#[non_exhaustive]
pub enum Error {
    /// A consumer depends on fields that the block it was given does not
    /// hold. A field is missing when nobody asked the reader for it, and
    /// when its source has none to give: a source built from an array of
    /// genotypes has no alleles. The consumer finds which ones with
    /// `asked_for.difference(block.fields())`, and one error names them
    /// all, so that a consumer that depends on two fields reports both.
    #[error("the block does not hold the fields this needs: {fields}")]
    FieldsNotInTheBlock {
        /// The fields that were asked for and that the block does not
        /// hold.
        fields: Needs,
    },

    /// A reader that takes a size was asked for blocks of 0 variants. A
    /// block holds one variant at least, and the caller that wants the
    /// size popnei chooses asks for none instead of asking for 0.
    #[error("a reader was asked for blocks of 0 variants, and a block holds 1 variant at least")]
    BlockOfNoVariants,

    /// A block of that many variants needs more memory than the machine
    /// gives: its genotypes, the variants times the individuals times the
    /// ploidy, are more than a `usize` holds, or one of its columns was
    /// asked of the machine and not given.
    ///
    /// The size is the one a caller asked for, or the one popnei chose for
    /// the individuals of the source when they asked for none: the VCF
    /// reader refuses the first when it is built and finds the second when
    /// it builds its first block, so that a file opened for its individuals
    /// alone is never refused for a size that nobody asked for. `size` says
    /// which of the two it is, and the message ends with what the caller
    /// does about it.
    #[error(
        "a block of {num_vars_per_block} variants of {num_individuals} individuals of the ploidy {ploidy} needs more memory than this machine gives; {way_out}",
        way_out = size.way_out()
    )]
    BlockTooLarge {
        /// How many variants a block was asked to hold.
        num_vars_per_block: usize,
        /// How many individuals the source has.
        num_individuals: usize,
        /// How many alleles the genotype of one individual holds.
        ploidy: usize,
        /// Whether that size is the one the caller asked for or the one
        /// popnei chose for these individuals.
        size: BlockSize,
    },

    /// A reader gave a block of other individuals or of another ploidy
    /// than it says its source has, so the rows of its blocks are not rows
    /// of one array of variants x individuals x ploidy and cannot be
    /// joined. `reblock` finds it, and it is a defect of that reader.
    #[error(
        "the reader says its source has {num_individuals} individuals of the ploidy {ploidy} and gave a block of {found_num_individuals} individuals of the ploidy {found_ploidy}"
    )]
    BlocksDoNotFitTogether {
        /// How many individuals the reader says its source has.
        num_individuals: usize,
        /// The ploidy the reader says its source has.
        ploidy: usize,
        /// How many individuals the block it gave has.
        found_num_individuals: usize,
        /// How many alleles the genotype of one individual holds in it.
        found_ploidy: usize,
    },

    /// A reader gave a block of no variants, which no reader of popnei
    /// does: every block a reader gives holds one variant at least, and a
    /// reader with no more variants gives no block. `reblock` finds it and
    /// gives no block after it, because a source that gives one says
    /// nothing about whether the variants that follow are there, and a
    /// reader that asked again would never come back from a source that
    /// always gives one.
    #[error(
        "a reader gave a block of no variants, and every block holds 1 variant at least; the reader that gave it has a defect"
    )]
    ReaderGaveABlockOfNoVariants,

    /// An array of a block is not of the size the block says: its
    /// genotypes are not its variants times its individuals times its
    /// ploidy, or a column has not one entry for each variant. The fields
    /// of a block are public, so a reader with a defect can build one, and
    /// its genotypes would be read one at the place of another with
    /// nothing to show it. `Block::check` is what finds it.
    #[error(
        "the `{array}` of a block holds {found} entries, and a block of its size holds {expected}"
    )]
    BlockArrayOfAnotherSize {
        /// The array that is not of the size of the block: `gts` or the
        /// name of a column.
        array: &'static str,
        /// How many entries it holds.
        found: usize,
        /// How many it has to hold.
        expected: usize,
    },

    /// A filter gave `Block::retain_vars` a number of values other than
    /// the variants of the block. There is one value for each variant, and
    /// the block is left as it was.
    #[error(
        "the variants to keep are {found} values and the block has {num_vars} variants; there is one value for each variant of the block"
    )]
    KeepOfAnotherSize {
        /// How many values were given.
        found: usize,
        /// How many variants the block holds.
        num_vars: usize,
    },

    /// A name that was given for a column of a block is not one of the
    /// five. It is a Python or a TypeScript user who writes them, in
    /// `iter_blocks(fields=...)`, so the message lists the names there are.
    #[error(
        "`{name}` is not a field of a block; the fields are {fields}",
        fields = crate::block::field_names_listed()
    )]
    NotAFieldOfABlock {
        /// The name that was given and is not a field of a block.
        name: String,
    },

    /// The source the VCF reader was given holds something else. A VCF
    /// starts with `#`, and a gzipped one with the two bytes of gzip.
    #[error("the source is not a VCF: it starts with {found}")]
    NotAVcf {
        /// The first bytes of the source, as text.
        found: String,
    },

    /// The header of the VCF is not one popnei can read. It needs the
    /// `#CHROM` line, its nine first columns and one individual or more
    /// after them, each with its own name.
    #[error("the header of the VCF cannot be read: {problem}")]
    VcfHeader {
        /// What is wrong with the header.
        problem: String,
    },

    /// The ploidy the VCF reader was asked for is 0, or above the largest
    /// one it reads. It is the one thing `VcfReader::new` refuses that does
    /// not come from the source.
    #[error(
        "the ploidy asked of the VCF reader is {ploidy}, and a genotype holds one allele at least and {largest} at most"
    )]
    VcfPloidyOutOfRange {
        /// The ploidy that was asked for.
        ploidy: usize,
        /// The largest one the reader takes, `vcf::MAX_PLOIDY`.
        largest: usize,
    },

    /// A data line of the VCF is not one popnei can read.
    #[error("line {line} of the VCF, {place}: {problem}")]
    VcfDataLine {
        /// The number of the line in the file, counted from 1 with the
        /// lines of the header.
        line: u64,
        /// The column or the individual the problem is in.
        place: VcfPlace,
        /// What is wrong there.
        problem: String,
    },

    /// A genotype of the VCF holds a number of alleles other than the
    /// ploidy the reader was given. popnei does not read a VCF of mixed
    /// ploidies: its calculations are not defined for one.
    #[error(
        "line {line} of the VCF, the column of {individual}: its genotype is of the ploidy {found} and the reader was asked for the ploidy {expected}; popnei does not read a VCF whose genotypes are of different ploidies, and the ploidy is an argument of the reader"
    )]
    VcfGenotypePloidy {
        /// The number of the line in the file, counted from 1 with the
        /// lines of the header.
        line: u64,
        /// The name of the individual whose genotype it is.
        individual: String,
        /// How many alleles the genotype holds.
        found: usize,
        /// The ploidy the reader was given.
        expected: usize,
    },

    /// The parse of a batch of lines of a VCF did not come back, which is
    /// what a panic inside it leaves behind: the lines of that batch were
    /// never parsed, and a reader that went on would give the variants of
    /// the lines that were and drop the others without a word. Nothing a
    /// VCF holds panics the parse, so this says that popnei has a defect,
    /// and the panic itself is what names it.
    #[error(
        "the parse of the lines of the VCF up to the line {line} did not come back, which a panic inside it leaves behind; the reader gives no more variants, because the lines it did not parse would be dropped without a word, and the panic that came before this error is what says where the defect is"
    )]
    VcfParseNotFinished {
        /// The number of the last line that was read, counted from 1 with
        /// the lines of the header.
        line: u64,
    },

    /// The source of the VCF was written by bgzip and does not end with the
    /// empty member of 28 bytes that marks the end of such a file, so its
    /// last bytes are missing: a download that stopped, a copy that was cut
    /// short. It comes where the reader would have said that there are no
    /// more variants.
    ///
    /// A gzip file that bgzip did not write has no such mark and is read to
    /// its end.
    #[error(
        "the VCF was written by bgzip and does not end with the empty member of 28 bytes that marks the end of a bgzipped file, so the file is cut short and the variants after the cut are not in it; the file has to be fetched or copied again. bcftools says of the same file `no BGZF EOF marker; file may be truncated`"
    )]
    VcfBgzipEndMissing,

    /// A member of the source of the VCF, which bgzip wrote, is not one
    /// bgzip could have written: its header is not that of a member of such
    /// a file, the size it states is not the size it has, the text that
    /// came out of it is not the text its CRC32 and its length describe, or
    /// the file goes on after the member that marks its end. The file was
    /// damaged after it was written, by a copy or a transfer that did not
    /// check what it carried.
    ///
    /// The member is counted from 1 and `offset` is the byte of the
    /// compressed file where it starts, so that a user can look at it with
    /// `xxd -s`.
    #[error(
        "the VCF was written by bgzip and its member {member}, which starts at the byte {offset} of the compressed file, is corrupted, so the file has to be fetched or copied again: {problem}"
    )]
    VcfBgzipCorrupted {
        /// Which member of the file it is, counted from 1.
        member: u64,
        /// The byte of the compressed file where that member starts,
        /// counted from 0.
        offset: u64,
        /// What is wrong with it.
        problem: String,
    },

    /// The file of a VCF, or of another source of variants, could not be
    /// opened. It carries the path, which `std::io::Error` does not, so
    /// that a message names the file and a binding can put it where its
    /// language keeps it, `OSError.filename` in Python.
    #[error("the file {path} could not be opened: {source}")]
    FileNotOpened {
        /// The path that was asked for.
        path: PathBuf,
        /// Why the file could not be opened.
        source: std::io::Error,
    },

    /// The bytes of a source could not be read.
    #[error("the source could not be read: {0}")]
    Io(#[from] std::io::Error),
}

/// What every operation of popnei that can fail returns.
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::Error;
    use crate::io::vcf::VcfPlace;
    use crate::variant::Needs;

    /// The message has to name the fields, because that is what tells the
    /// caller which reader to ask or which calculation to drop.
    #[test]
    fn the_message_of_a_field_the_block_does_not_hold_names_the_field() {
        let error = Error::FieldsNotInTheBlock {
            fields: Needs::ALLELES | Needs::QUAL,
        };
        let message = error.to_string();
        assert!(message.contains("alleles"), "{message}");
        assert!(message.contains("qual"), "{message}");
        assert!(!message.contains("gts"), "{message}");
    }

    /// A user who gets one of these has the file open in front of them, so
    /// the message says which line and which column or individual to look
    /// at.
    #[test]
    fn the_message_of_a_wrong_data_line_names_the_line_and_the_place() {
        let error = Error::VcfDataLine {
            line: 12,
            place: VcfPlace::Column("POS"),
            problem: "`x` is not a position".to_string(),
        };
        let message = error.to_string();
        assert!(message.contains("12"), "{message}");
        assert!(message.contains("POS"), "{message}");
        assert!(message.contains("`x` is not a position"), "{message}");

        let error = Error::VcfGenotypePloidy {
            line: 9,
            individual: "ind2".to_string(),
            found: 4,
            expected: 2,
        };
        let message = error.to_string();
        assert!(message.contains("line 9"), "{message}");
        assert!(message.contains("ind2"), "{message}");
        assert!(message.contains("ploidy 4"), "{message}");
        assert!(message.contains("ploidy 2"), "{message}");
        // The reader can be asked for another ploidy, and a VCF of mixed
        // ploidies is refused whatever it is asked for: a user who gets
        // this needs to be told both.
        assert!(message.contains("different ploidies"), "{message}");
        assert!(message.contains("argument"), "{message}");
    }
}
