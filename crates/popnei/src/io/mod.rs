//! The files the variants of popnei come from and go to.
//!
//! [`vcf`] reads a VCF, plain or gzipped, and gives its variants in blocks,
//! through the [`BlockReader`](crate::block::BlockReader) trait, which is
//! how every source of variants of popnei gives them. `bgzf`, which is of
//! the crate and not of its public documentation, is what it reads a source
//! that bgzip wrote with, member by member. The reader and
//! the writer of the vars file, popnei's own file of variants, are not
//! written yet; section 6 of `docs/architecture.md` has its format.

pub(crate) mod bgzf;
pub mod vcf;
