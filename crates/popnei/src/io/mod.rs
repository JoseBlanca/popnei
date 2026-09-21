//! The files the variants of popnei come from and go to.
//!
//! [`vcf`] reads a VCF, plain or gzipped, and gives its variants in blocks,
//! through the [`BlockReader`](crate::block::BlockReader) trait, which is
//! how every source of variants of popnei gives them. `bgzf`, which is of
//! the crate and not of its public documentation, is what it reads a source
//! that bgzip wrote with, member by member.
//!
//! [`vars`] is the vars file, popnei's own file of variants: one arrow IPC
//! file with a record batch for each block, which any program with an arrow
//! library opens as a table. What is there of it is what the file says about
//! itself, the two keys whose values are json; its writer and its reader are
//! being written. `docs/specs/io_vars.md` has the format and section 6 of
//! `docs/architecture.md` where it sits.

pub(crate) mod bgzf;
pub mod vars;
pub mod vcf;
