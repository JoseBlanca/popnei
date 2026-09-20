//! The files the variants of popnei come from and go to.
//!
//! [`vcf`] reads a VCF, plain or gzipped, and gives its variants in blocks,
//! through the [`BlockReader`](crate::block::BlockReader) trait. The reader
//! and the writer of the vars file, popnei's own file of variants, are not
//! written yet; section 6 of `docs/architecture.md` has its format.

pub mod vcf;
