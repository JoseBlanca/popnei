//! What the two sources of variants share: one pass over a source, the
//! columns of the blocks that pass gives JavaScript, and the bytes of the
//! vars file written from it.
//!
//! A source is the bytes of a file with what is needed to read it, a VCF
//! with its options in `vcf.rs` and a vars file in `vars.rs`. Each is a
//! class of wasm-bindgen that the `Variants` of the TypeScript package
//! holds, each reads what says what the file holds when it is built, so
//! bytes that are not of its format fail at `openVcf` or at `openVars`, and
//! each reads the bytes again at every pass, which is what lets a user give
//! the same `Variants` to one calculation after another. What they have in
//! common is [`OpenSource`], the reader of one pass.
//!
//! [`Blocks`] is that pass, whichever source it came from: it owns a reader
//! of blocks of the core with a `Reblock` at its end, which gives the blocks
//! the size that was asked for. [`BlockColumns`] is one block on its way out.
//!
//! The three live in the memory of wasm, which the garbage collector of
//! JavaScript does not see, so the TypeScript package frees each of them:
//! the block as soon as its columns are read, the pass when the iteration
//! ends or is abandoned, and the source when the user calls `free()`.
//!
//! The columns leave as JavaScript arrays, which wasm-bindgen copies out of
//! the memory of wasm: a typed array that is a view into that memory stops
//! being valid when it grows, as section 11 of `docs/architecture.md` says.
//! Each column is moved out of the block as it is read, so that the copy
//! that crosses is the only one.

use std::io::{Cursor, ErrorKind, Write};
use std::sync::Arc;

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::{AllelesColumn, Block, BlockReader, Reblock, needs_of_the_fields};
use popnei::variant::Needs;

use crate::errors::JsPopneiError;

/// The largest position a block hands to JavaScript, 2^53.
///
/// The positions cross as float64, which holds every whole number up to
/// this one and not the ones above it: 2^53 + 1 would arrive as 2^53.
const LARGEST_POSITION: u64 = 9_007_199_254_740_992;

/// A file of variants that was opened, which every pass reads again.
pub(crate) trait OpenSource {
    /// The reader of one pass over the source, which reads the bytes from
    /// their start.
    ///
    /// `num_vars_per_block` is the size the caller will ask the blocks for,
    /// which a source whose reader can give them at that size is built
    /// with, so that the `Reblock` over it has nothing to cut or to join. A
    /// source that gives its blocks as they are in the file, the vars file
    /// whose batches were written at one size, leaves it to that `Reblock`.
    ///
    /// # Errors
    ///
    /// When what says what the file holds, the header of a VCF or the
    /// schema of a vars file, cannot be read.
    fn reader(
        &self,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error>;
}

/// The bytes of a file, which every pass over it shares.
///
/// `Cursor` reads anything that gives a `&[u8]`, and an `Arc<Vec<u8>>` gives
/// a `&Vec<u8>`, so this is the one line that turns the second into the
/// first. What it holds is the `Vec` wasm-bindgen filled with the bytes of
/// the `Uint8Array`, and nothing copies it again: a copy of a file of 80 MB
/// stays in the memory of wasm for as long as the tab lives, because that
/// memory grows and never shrinks.
pub(crate) struct SharedBytes(Arc<Vec<u8>>);

impl AsRef<[u8]> for SharedBytes {
    fn as_ref(&self) -> &[u8] {
        self.0.as_slice()
    }
}

/// One pass over `bytes`, from their first byte, which shares them with
/// every other pass over the same source.
pub(crate) fn cursor_of(bytes: &Arc<Vec<u8>>) -> Cursor<SharedBytes> {
    Cursor::new(SharedBytes(Arc::clone(bytes)))
}

/// One pass over `source`, whose blocks hold `fields` besides the genotypes,
/// `num_vars_per_block` variants each.
///
/// # Errors
///
/// When a name of `fields` is not a field of a block, when
/// `num_vars_per_block` is 0, when the source cannot be read again, and when
/// the genotypes of one block are more than the memory of wasm addresses.
pub(crate) fn blocks_of(
    source: &dyn OpenSource,
    fields: Vec<String>,
    num_vars_per_block: Option<usize>,
) -> Result<Blocks, JsPopneiError> {
    let needs = needs_of_the_fields(fields.iter().map(String::as_str))?;
    // The source is asked for the size the user wants, so a reader that can
    // give it has nothing for the `Reblock` over it to cut or to join and
    // every block goes through with no copy. That `Reblock` is there for the
    // sources that give another size, the vars file whose batches were
    // written at one size and a filter among them, and it is what
    // `docs/specs/block.md` puts at the end of every `iterBlocks`.
    let mut reader = source.reader(num_vars_per_block)?;
    reader.set_needs(needs.union(Needs::GTS));
    Ok(Blocks {
        reader: Box::new(Reblock::new(reader, num_vars_per_block)?),
        finished: false,
    })
}

/// How many bytes one piece of a vars file that is being written holds,
/// 1 MiB.
///
/// It is what the memory of wasm grows by at a time while a file is written
/// and what one call of `next_piece` copies out, and nothing of the format
/// depends on it: the pieces are put together in JavaScript into the one
/// array the user gets. Writing a file of 18.3 MB from a vars source grew
/// that memory by 33.3 MB with pieces of this size, by 34.2 MB with pieces
/// of 128 KiB and by 36.8 MB with pieces of 4 MiB, so what is left above
/// the file is not the pieces but what the reader and arrow-rs hold while a
/// batch is written.
const BYTES_PER_PIECE: usize = 1024 * 1024;

/// A vars file that was written, held in pieces of [`BYTES_PER_PIECE`].
///
/// The whole file is in the memory of wasm when the write is over, and that
/// memory never shrinks, so what a `Vec` that grew by doubling cost a tab
/// was the bytes it had and the buffer twice that size it copied them into,
/// at every doubling, all of it kept for as long as the page lives.
/// Measured on a vars file of 18.3 MB written from a vars source, that one
/// grew the memory of wasm by 62.4 MB and these pieces grow it by 33.3 MB.
///
/// Every piece leaves the memory of wasm as it is read, which is where the
/// copy that crosses is made, and the package puts them together into the
/// `Uint8Array` the user gets, in the heap of JavaScript.
#[wasm_bindgen]
pub struct VarsFile {
    /// The pieces in the order they were written, each of them empty once
    /// it has been given to JavaScript.
    pieces: Vec<Vec<u8>>,
    num_bytes: usize,
    /// Which piece is the next one to give.
    next: usize,
}

#[wasm_bindgen]
impl VarsFile {
    /// How many bytes the whole file holds.
    #[must_use]
    pub fn num_bytes(&self) -> usize {
        self.num_bytes
    }

    /// The next piece of the file, or `undefined` when it has all been
    /// given. Each of them leaves the memory of wasm as it is read.
    pub fn next_piece(&mut self) -> Option<Vec<u8>> {
        let piece = self.pieces.get_mut(self.next)?;
        self.next = self.next.saturating_add(1);
        Some(std::mem::take(piece))
    }
}

/// The sink the core writes a vars file into: it takes the bytes in pieces
/// of [`BYTES_PER_PIECE`] and never copies what it has into a larger buffer.
struct PiecesOfTheFile {
    pieces: Vec<Vec<u8>>,
    num_bytes: usize,
}

impl PiecesOfTheFile {
    fn new() -> PiecesOfTheFile {
        PiecesOfTheFile {
            pieces: Vec::new(),
            num_bytes: 0,
        }
    }

    /// The piece the next bytes go into, a new one when the last is full.
    ///
    /// # Errors
    ///
    /// When the memory of wasm cannot take another piece, which is the
    /// error of a vars file that could not be written: the core wraps what
    /// a sink says, and a tab that has no memory left for the file is told
    /// so instead of trapping on a failed allocation.
    fn piece_with_room(&mut self) -> std::io::Result<&mut Vec<u8>> {
        let full = match self.pieces.last() {
            Some(piece) => piece.len() >= piece.capacity(),
            None => true,
        };
        if full {
            let mut piece: Vec<u8> = Vec::new();
            piece.try_reserve_exact(BYTES_PER_PIECE).map_err(|_| {
                std::io::Error::new(
                    ErrorKind::OutOfMemory,
                    format!(
                        "the memory of this tab does not take {BYTES_PER_PIECE} bytes \
                         more of the file, which holds {num_bytes} bytes so far",
                        num_bytes = self.num_bytes
                    ),
                )
            })?;
            self.pieces.push(piece);
        }
        self.pieces
            .last_mut()
            .ok_or_else(|| std::io::Error::other("the file has no piece to be written into"))
    }
}

impl Write for PiecesOfTheFile {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.is_empty() {
            return Ok(0);
        }
        let piece = self.piece_with_room()?;
        let room = piece.capacity().saturating_sub(piece.len());
        let taken = room.min(bytes.len());
        let Some(head) = bytes.get(..taken) else {
            return Err(std::io::Error::other(
                "the bytes of the file are fewer than what is being taken from them",
            ));
        };
        piece.extend_from_slice(head);
        self.num_bytes = self.num_bytes.checked_add(taken).ok_or_else(|| {
            std::io::Error::other("the file holds more bytes than this machine counts")
        })?;
        Ok(taken)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Every variant of `source` as a vars file, one batch of
/// `num_vars_per_block` variants after another, and `None` for the size
/// popnei chooses for the individuals of the source.
///
/// The file is built in the memory of wasm, in pieces that cross one by one:
/// a `Uint8Array` that were a view into that memory would stop being valid
/// the next time it grows. A tab holds the source and the file at once, so
/// the memory it needs is the two together.
///
/// # Errors
///
/// When `num_vars_per_block` is 0, when the source cannot be read, when a
/// block of it is not one a vars file holds, and when the memory of the tab
/// does not take the file.
pub(crate) fn bytes_of_a_vars_file(
    source: &dyn OpenSource,
    num_vars_per_block: Option<usize>,
) -> Result<VarsFile, JsPopneiError> {
    // The source is asked for the size the batches will have, as a pass is,
    // so that the `reblock` the core puts over it has nothing to cut or to
    // join. What that saves is the memory of a block: a VCF read with the
    // size popnei chooses for 1000 individuals, 5000 variants, holds 10.4 MB
    // of genotypes while it is written, and 0.3 MB when the caller asked for
    // batches of 100. A source that cannot give that size, the vars file
    // whose batches were written at another one, leaves it to the `reblock`.
    let reader = source.reader(num_vars_per_block)?;
    let written = popnei::io::vars::write_vars(reader, PiecesOfTheFile::new(), num_vars_per_block)?;
    Ok(VarsFile {
        pieces: written.pieces,
        num_bytes: written.num_bytes,
        next: 0,
    })
}

/// One pass over a source of variants, which gives them block by block.
#[wasm_bindgen]
pub struct Blocks {
    reader: Box<dyn BlockReader>,
    /// Whether the pass is over: the reader has no more blocks, or a block
    /// was lost with an error. After either there is no block.
    finished: bool,
}

#[wasm_bindgen]
impl Blocks {
    /// The next block of the pass, or `undefined` when the source has no
    /// more variants.
    ///
    /// The names of the chromosomes are taken after the block was given, as
    /// `docs/specs/block.md` says: the table of the reader grows while the
    /// file is read, and a block holds numbers of it.
    ///
    /// The block is checked before its genotypes cross as one `Int8Array`
    /// of variants x individuals x ploidy: a block whose arrays are not of
    /// its size would be read one genotype at the place of another, with
    /// nothing to show it.
    ///
    /// # Errors
    ///
    /// When a variant cannot be read, when the block is not of its own size,
    /// and when a position of the block is above [`LARGEST_POSITION`]. The
    /// block that was being built is lost with the error, and every call
    /// after it gives no block: the errors of this pass happen after the
    /// block was taken from the reader, which knows nothing of them, so it
    /// is this pass that keeps the promise of `docs/specs/block.md` that
    /// the variants after a wrong one are not handed out as if nothing had
    /// happened.
    pub fn next_block(&mut self) -> Result<Option<BlockColumns>, JsPopneiError> {
        if self.finished {
            return Ok(None);
        }
        let columns = self.columns_of_the_next_block();
        if !matches!(columns, Ok(Some(_))) {
            self.finished = true;
        }
        columns
    }
}

impl Blocks {
    /// The columns of the next block of the reader, or `None` when it has no
    /// more variants. What ends the pass is [`Blocks::next_block`], which
    /// calls this one.
    fn columns_of_the_next_block(&mut self) -> Result<Option<BlockColumns>, JsPopneiError> {
        let Some(block) = self.reader.next_block()? else {
            return Ok(None);
        };
        block.check()?;
        let chroms = self.reader.chroms();
        let names = |numbers: Vec<u32>| {
            numbers
                .into_iter()
                .map(|number| {
                    chroms.name(number).map(str::to_owned).ok_or_else(|| {
                        JsPopneiError::Broken(format!(
                            "the chromosome number {number} of a block is not in the \
                             table of the reader that gave it"
                        ))
                    })
                })
                .collect::<Result<Vec<String>, JsPopneiError>>()
        };
        let Block {
            num_vars,
            num_individuals,
            ploidy,
            gts,
            chrom,
            pos,
            id,
            alleles,
            qual,
        } = block;
        let alleles = alleles.map(|column| alleles_of(&column)).transpose()?;
        let (alleles, num_alleles_per_var) = match alleles {
            Some((texts, counts)) => (Some(texts), Some(counts)),
            None => (None, None),
        };
        Ok(Some(BlockColumns {
            num_vars,
            num_individuals,
            ploidy,
            gts: Some(gts),
            chrom: chrom.map(names).transpose()?,
            pos: pos.map(positions_of).transpose()?,
            id,
            alleles,
            num_alleles_per_var,
            qual,
        }))
    }
}

/// The columns of one block on their way to TypeScript.
///
/// Each of them leaves the block the first time it is asked for, and the
/// call after that gives nothing: the package reads every column once, into
/// the object a user holds, and frees the block. What a column costs is then
/// one copy, the one wasm-bindgen makes when it crosses.
#[wasm_bindgen]
pub struct BlockColumns {
    num_vars: usize,
    num_individuals: usize,
    ploidy: usize,
    gts: Option<Vec<i8>>,
    chrom: Option<Vec<String>>,
    pos: Option<Vec<f64>>,
    id: Option<Vec<String>>,
    /// The alleles of every variant one variant after another, and how many
    /// of them each variant has: a JavaScript array of arrays of strings is
    /// not one of the types wasm-bindgen carries, so the column crosses as
    /// the core holds it, flat, and the package cuts it into the alleles of
    /// each variant.
    alleles: Option<Vec<String>>,
    num_alleles_per_var: Option<Vec<u32>>,
    qual: Option<Vec<f32>>,
}

#[wasm_bindgen]
impl BlockColumns {
    /// How many variants the block holds.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals the source has, the same for every variant.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many alleles the genotype of one individual holds.
    #[must_use]
    pub fn ploidy(&self) -> usize {
        self.ploidy
    }

    /// The genotypes, `num_vars` x `num_individuals` x `ploidy` alleles,
    /// variant after variant and inside a variant individual after
    /// individual. Every block holds them, so `undefined` says that they
    /// were read already, as it does for every other column.
    pub fn gts(&mut self) -> Option<Vec<i8>> {
        self.gts.take()
    }

    /// The name of the chromosome of each variant, or `undefined` when the
    /// pass was not asked for it.
    pub fn chrom(&mut self) -> Option<Vec<String>> {
        self.chrom.take()
    }

    /// The position of each variant, 1 based as in a VCF.
    pub fn pos(&mut self) -> Option<Vec<f64>> {
        self.pos.take()
    }

    /// The id of each variant, an empty text for a variant that has none.
    pub fn id(&mut self) -> Option<Vec<String>> {
        self.id.take()
    }

    /// The alleles of every variant, the reference one first, one variant
    /// after another: the first `num_alleles_per_var()[0]` of them are the
    /// alleles of the first variant, and so on.
    pub fn alleles(&mut self) -> Option<Vec<String>> {
        self.alleles.take()
    }

    /// How many alleles each variant has, the reference one counted.
    pub fn num_alleles_per_var(&mut self) -> Option<Vec<u32>> {
        self.num_alleles_per_var.take()
    }

    /// The quality of each variant, phred scaled as the QUAL of a VCF, and
    /// NaN for a variant that has none.
    pub fn qual(&mut self) -> Option<Vec<f32>> {
        self.qual.take()
    }
}

/// The positions of a block as the numbers of JavaScript, which are float64
/// and hold a position of up to 2^53 exactly, as `docs/specs/block.md` says.
///
/// # Errors
///
/// When a position is above [`LARGEST_POSITION`]. A float64 would round it,
/// and the same file read from Python gives the position the source has, so
/// the two languages would disagree about where a variant is. No genome
/// comes near that number: the longest chromosome that has been assembled
/// is 2.5e8 bases.
fn positions_of(positions: Vec<u64>) -> Result<Vec<f64>, JsPopneiError> {
    positions
        .into_iter()
        .map(|pos| {
            if pos > LARGEST_POSITION {
                return Err(JsPopneiError::NotInJavaScript(format!(
                    "the position {pos} of a variant is above {LARGEST_POSITION}, the \
                     largest whole number a number of JavaScript holds"
                )));
            }
            Ok(pos as f64)
        })
        .collect()
}

/// The alleles of every variant of a block, the reference one first, one
/// variant after another, and how many of them each variant has.
///
/// # Errors
///
/// When a variant has more alleles than a `u32` counts, which no source
/// gives: a variant of a VCF has as many alleles as its ALT column names.
fn alleles_of(column: &AllelesColumn) -> Result<(Vec<String>, Vec<u32>), JsPopneiError> {
    let mut texts = Vec::new();
    let mut counts = Vec::with_capacity(column.num_vars());
    for var in 0..column.num_vars() {
        let num_alleles = column.num_alleles(var);
        texts.extend((0..num_alleles).map(|allele| column.allele(var, allele).to_owned()));
        counts.push(u32::try_from(num_alleles).map_err(|_| {
            JsPopneiError::Broken(format!(
                "the variant {var} of a block has {num_alleles} alleles, more than a \
                 JavaScript array of counts holds"
            ))
        })?);
    }
    Ok((texts, counts))
}
