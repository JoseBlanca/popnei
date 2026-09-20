//! The error of the crate.
//!
//! Every module of popnei adds its cases to [`Error`], and every operation
//! that can fail returns the [`Result`] of this module. The enum is
//! `#[non_exhaustive]`, so a module that is written later adds a case
//! without breaking the code that matches on it.

use thiserror::Error as ThisError;

use crate::io::vcf::VcfPlace;
use crate::variant::Needs;

/// Anything that went wrong in popnei.
#[derive(Debug, ThisError)]
#[non_exhaustive]
pub enum Error {
    /// The consumer depends on fields that the reader did not fill. A
    /// reader may leave out a field that was asked for when its source has
    /// none, an array of genotypes that has no alleles, and the consumer
    /// finds it in `filled` of the variant.
    #[error("the reader did not fill the fields that were asked for: {fields}")]
    FieldsNotFilled {
        /// The fields that were asked for and are not in `filled`.
        fields: Needs,
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
    fn the_message_of_a_field_that_was_not_filled_names_the_field() {
        let error = Error::FieldsNotFilled {
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
