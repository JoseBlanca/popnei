//! The error of the crate.
//!
//! Every module of popnei adds its cases to [`Error`], and every operation
//! that can fail returns the [`Result`] of this module. The enum is
//! `#[non_exhaustive]`, so a module that is written later adds a case
//! without breaking the code that matches on it.

use thiserror::Error as ThisError;

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
}

/// What every operation of popnei that can fail returns.
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::Error;
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
}
