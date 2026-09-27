//! The members of bgzip that a bgzipped VCF is written in, as "How it runs"
//! of the VCF writer in `docs/specs/io_vcf.md` has them: members of 65280
//! bytes of text, each a gzip stream whose header states its size in the
//! extra field `BC`, compressed on the threads of rayon and written in the
//! order of the text, and the empty member that marks the end of the file.

use std::io::Write;

use flate2::{Compress, CompressError, Compression, Crc, FlushCompress, Status};

use super::not_written;
use crate::error::{Error, Result};

/// How many bytes of text a member holds, the last one of a file aside:
/// 65280, what bgzip fills one with, which leaves room in the 65536 bytes a
/// member has for a text that deflate cannot shrink.
pub(crate) const TEXT_OF_A_MEMBER: usize = 65280;

/// The most bytes of a member, header and end included: its size is stated
/// in two bytes as the size less 1.
const MOST_BYTES_OF_A_MEMBER: usize = 65536;

/// The header of a member before the two bytes of its size: the two bytes
/// of gzip, the method deflate, the flag of an extra field, no time stamp,
/// no extra flags, the system 255, unknown, as bgzip writes it, the six
/// bytes of the extra field, and the name and the length of its subfield
/// `BC`.
const HEADER_BEFORE_THE_SIZE: [u8; 16] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
];

/// The bytes of the header, the two of the size included, and the bytes
/// after the data: the CRC32 of the text and its length.
const BYTES_OF_THE_HEADER: usize = 18;
const BYTES_AFTER_THE_DATA: usize = 8;

/// The most bytes of compressed data a member holds.
const MOST_DATA_OF_A_MEMBER: usize =
    MOST_BYTES_OF_A_MEMBER - BYTES_OF_THE_HEADER - BYTES_AFTER_THE_DATA;

/// The empty member that ends a file bgzip wrote, the 28 bytes htslib
/// writes: a member of no text, whose deflate is the two bytes `03 00`.
pub(crate) const THE_EMPTY_MEMBER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// The level of deflate of every member: 6, zlib's default, which bgzip
/// 1.24 uses when it is given none. On `big.vcf`, 100000 variants x 1000
/// individuals, on one thread of the owner's M5 Pro on 27 September 2026,
/// `write_vcf` bgzipped with flate2 over zlib-rs took 5.48 s at this level
/// and gave 37.4 MB, and `bcftools view -Oz` 8.84 s and 37.7 MB. With
/// miniz_oxide, the backend before, it took 10.66 s and gave 38.4 MB, and
/// on 26 September zlib-rs at level 5 took 3.1 s and gave 39.8 MB. The
/// owner chose zlib-rs at 6 on 27 September 2026.
pub(crate) const COMPRESSION_LEVEL: u32 = 6;

/// The first byte of a deflate block that stores its text as it is and is
/// the last block of its stream.
const THE_LAST_STORED_BLOCK: u8 = 0x01;

/// The text of a bgzipped file on its way to the sink: the members are cut
/// straight out of the text of each block, where the writer formatted it,
/// and what does not fill a member waits here for the text of the next.
pub(super) struct BgzipOut<W: Write> {
    sink: W,
    /// The text that did not fill a member, less than 65280 bytes.
    waiting: Vec<u8>,
    /// One buffer for each member compressed at once, kept from one block
    /// to the next.
    members: Vec<Vec<u8>>,
}

/// The pieces of text one member holds, in their order, 65280 bytes in all
/// but for the last member of the file: the tail of the text that waited,
/// and the parts of the buffers of a block that follow it.
type TextOfAMember<'a> = Vec<&'a [u8]>;

impl<W: Write> BgzipOut<W> {
    pub(super) fn new(sink: W) -> BgzipOut<W> {
        BgzipOut {
            sink,
            waiting: Vec::new(),
            members: Vec::new(),
        }
    }

    /// The text of `pieces`, in their order, after the text given before:
    /// every member it fills is compressed and written, and what fills none
    /// waits.
    ///
    /// # Errors
    ///
    /// When the sink refuses a member, and when a member could not be put
    /// together, which is a defect.
    pub(super) fn write_the_text(&mut self, pieces: &[&[u8]]) -> Result<()> {
        let BgzipOut {
            sink,
            waiting,
            members,
        } = self;
        let (full, tail) = members_of(waiting, pieces);
        if !full.is_empty() {
            members.resize_with(full.len(), Vec::new);
            members.truncate(full.len());
            compress_the_members(&full, members)?;
            for member in members.iter() {
                sink.write_all(member).map_err(not_written)?;
            }
        }
        let mut left: Vec<u8> = Vec::new();
        for piece in &tail {
            left.extend_from_slice(piece);
        }
        *waiting = left;
        Ok(())
    }

    /// The member of the text that is left, when there is any, the empty
    /// member of the end, and the sink back with every byte flushed to it.
    ///
    /// # Errors
    ///
    /// When the sink refuses them.
    pub(super) fn finish(mut self) -> Result<W> {
        if !self.waiting.is_empty() {
            let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
            let mut member = Vec::new();
            compress_a_member(&mut compress, &[&self.waiting], &mut member)?;
            self.sink.write_all(&member).map_err(not_written)?;
        }
        self.sink
            .write_all(&THE_EMPTY_MEMBER)
            .map_err(not_written)?;
        self.sink.flush().map_err(not_written)?;
        Ok(self.sink)
    }
}

/// The text of `waiting` and then of `pieces`, cut into the members of
/// 65280 bytes it fills, each the pieces of it that it holds, and the
/// pieces of what is left, which fill no member. Nothing is copied.
fn members_of<'a>(
    waiting: &'a [u8],
    pieces: &[&'a [u8]],
) -> (Vec<TextOfAMember<'a>>, TextOfAMember<'a>) {
    let mut full = Vec::new();
    let mut current: TextOfAMember<'a> = Vec::new();
    let mut in_the_current = 0usize;
    for piece in std::iter::once(waiting).chain(pieces.iter().copied()) {
        let mut rest = piece;
        while !rest.is_empty() {
            let room = TEXT_OF_A_MEMBER.saturating_sub(in_the_current);
            let (taken, after) = rest.split_at(rest.len().min(room));
            current.push(taken);
            in_the_current = in_the_current.saturating_add(taken.len());
            rest = after;
            if in_the_current == TEXT_OF_A_MEMBER {
                full.push(std::mem::take(&mut current));
                in_the_current = 0;
            }
        }
    }
    (full, current)
}

/// The member of each of `texts` into `members`, one for each, in their
/// order. Natively they are compressed on the threads of rayon, with one
/// deflate for each job rayon splits them into, which on 18 threads is
/// about one for each member: its memory is a few hundred KB, which is
/// small beside the text of a member, and a deflate for each thread has
/// not been timed against it.
///
/// # Errors
///
/// What [`compress_a_member`] refuses.
#[cfg(not(target_family = "wasm"))]
fn compress_the_members(texts: &[TextOfAMember<'_>], members: &mut [Vec<u8>]) -> Result<()> {
    use rayon::iter::{
        IndexedParallelIterator, IntoParallelRefIterator, IntoParallelRefMutIterator,
        ParallelIterator,
    };

    texts
        .par_iter()
        .zip(members.par_iter_mut())
        .try_for_each_init(
            || Compress::new(Compression::new(COMPRESSION_LEVEL), false),
            |compress, (text, member)| compress_a_member(compress, text, member),
        )
}

/// The same members compressed one after another, which is what wasm does:
/// it has no threads.
///
/// # Errors
///
/// What [`compress_a_member`] refuses.
#[cfg(target_family = "wasm")]
fn compress_the_members(texts: &[TextOfAMember<'_>], members: &mut [Vec<u8>]) -> Result<()> {
    let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
    for (text, member) in texts.iter().zip(members.iter_mut()) {
        compress_a_member(&mut compress, text, member)?;
    }
    Ok(())
}

/// The member of the text of `pieces`, 65280 bytes at most, into `member`,
/// over what it held: the header with the size of the member, the text
/// deflated, and its CRC32 and length.
///
/// The text of several pieces is joined into one buffer before it is
/// deflated, since zlib-rs gives other bytes for a text it is fed in pieces
/// than for the same text fed whole, and the pieces are where the lines of
/// a block were formatted, so without the copy the file would depend on the
/// size of the blocks of its source. The copy is of 65280 bytes at most, of
/// about one member in four of `big.vcf`.
///
/// A deflate whose data do not fit the member, or whose stream did not end,
/// gives a member in which the writer stores the text as it is, one stored
/// block of deflate, which always fits and decompresses to the same text.
/// zlib-rs stores a text it cannot shrink on its own, 65326 bytes
/// of member for 65280 bytes of random text, so this is for a deflate that
/// does otherwise, which none of the tests reaches: [`the_deflate_fits`] is
/// tested on its own.
///
/// # Errors
///
/// When the text is longer than a member holds, which is a defect of the
/// caller.
fn compress_a_member(
    compress: &mut Compress,
    pieces: &[&[u8]],
    member: &mut Vec<u8>,
) -> Result<()> {
    let num_bytes = pieces
        .iter()
        .fold(0usize, |bytes, piece| bytes.saturating_add(piece.len()));
    let length = u16::try_from(num_bytes)
        .ok()
        .filter(|_| num_bytes <= TEXT_OF_A_MEMBER)
        .ok_or(Error::VcfWriterMemberNotBuilt {
            what: "text",
            found: num_bytes,
            most: TEXT_OF_A_MEMBER,
        })?;
    member.clear();
    member.reserve(MOST_BYTES_OF_A_MEMBER);
    member.extend_from_slice(&HEADER_BEFORE_THE_SIZE);
    member.extend_from_slice(&[0, 0]);
    let joined: Vec<u8>;
    let text: &[u8] = match pieces {
        [] => &[],
        [whole] => whole,
        _ => {
            joined = pieces.concat();
            &joined
        }
    };
    compress.reset();
    let deflated = deflate_the_text(compress, text, member);
    let data = member.len().saturating_sub(BYTES_OF_THE_HEADER);
    if !the_deflate_fits(&deflated, data) {
        member.truncate(BYTES_OF_THE_HEADER);
        store(pieces, length, member);
    }
    end_the_member(pieces, length, member)
}

/// The text deflated into `member` as one stream, and the status the
/// deflate gave at its end. A text that the deflate did not take whole,
/// which is what it does when the room of the member runs out, gives the
/// status `BufError`, and the writer stores the text then.
fn deflate_the_text(
    compress: &mut Compress,
    text: &[u8],
    member: &mut Vec<u8>,
) -> std::result::Result<Status, CompressError> {
    let status = compress.compress_vec(text, member, FlushCompress::Finish);
    match u64::try_from(text.len()).ok() == Some(compress.total_in()) {
        true => status,
        false => status.and(Ok(Status::BufError)),
    }
}

/// Whether a deflate that gave `deflated` and `data` bytes goes into the
/// member: its stream ended, and its data leave room in the 65536 bytes of
/// a member for the header and the end, 65510 bytes at most. Otherwise the
/// writer stores the text itself.
fn the_deflate_fits(deflated: &std::result::Result<Status, CompressError>, data: usize) -> bool {
    matches!(deflated, Ok(Status::StreamEnd)) && data <= MOST_DATA_OF_A_MEMBER
}

/// The data of a member that stores the text of `pieces`, of `length` bytes,
/// as it is: one deflate block of no compression, after the header in
/// `member`.
fn store(pieces: &[&[u8]], length: u16, member: &mut Vec<u8>) {
    member.push(THE_LAST_STORED_BLOCK);
    member.extend_from_slice(&length.to_le_bytes());
    member.extend_from_slice(&(!length).to_le_bytes());
    for piece in pieces {
        member.extend_from_slice(piece);
    }
}

/// The CRC32 and the length of the text of `pieces` after the data in
/// `member`, and the size of the member in its header.
///
/// # Errors
///
/// When the member is more bytes than its size states in two bytes, which
/// the callers make impossible.
fn end_the_member(pieces: &[&[u8]], length: u16, member: &mut Vec<u8>) -> Result<()> {
    let mut crc = Crc::new();
    for piece in pieces {
        crc.update(piece);
    }
    member.extend_from_slice(&crc.sum().to_le_bytes());
    member.extend_from_slice(&u32::from(length).to_le_bytes());
    // The stored member is 65280 bytes of text and 31 of the rest at most,
    // and a deflated one fits by the test above, so its size less 1 fits in
    // two bytes.
    let size = member
        .len()
        .checked_sub(1)
        .and_then(|size| u16::try_from(size).ok())
        .ok_or(Error::VcfWriterMemberNotBuilt {
            what: "whole member",
            found: member.len(),
            most: MOST_BYTES_OF_A_MEMBER,
        })?;
    if let Some(slot) = member.get_mut(16..BYTES_OF_THE_HEADER) {
        slot.copy_from_slice(&size.to_le_bytes());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(
        clippy::arithmetic_side_effects,
        reason = "the sizes and the places of the members of the small files of the tests"
    )]

    use std::io::Read;

    use flate2::read::MultiGzDecoder;
    use flate2::{Compress, Compression};

    use flate2::Status;

    use crate::error::Error;

    use super::{
        BYTES_OF_THE_HEADER, COMPRESSION_LEVEL, HEADER_BEFORE_THE_SIZE, MOST_BYTES_OF_A_MEMBER,
        TEXT_OF_A_MEMBER, compress_a_member, end_the_member, members_of, store, the_deflate_fits,
    };

    /// The member of `text`, and the text it decompresses to.
    fn member_and_text(text: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
        let mut member = Vec::new();
        compress_a_member(&mut compress, &[text], &mut member).expect("the member");
        let mut back = Vec::new();
        MultiGzDecoder::new(&member[..])
            .read_to_end(&mut back)
            .expect("the member decompresses");
        (member, back)
    }

    /// 65280 bytes of a generator of xorshift, which deflate does not
    /// shrink.
    fn bytes_that_do_not_shrink() -> Vec<u8> {
        let mut state: u32 = 2463534242;
        (0..TEXT_OF_A_MEMBER)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                state.to_le_bytes()[0]
            })
            .collect()
    }

    #[test]
    fn a_member_of_text_zlib_rs_cannot_shrink_is_stored_by_zlib_rs_and_fits() {
        let text = bytes_that_do_not_shrink();
        let (member, back) = member_and_text(&text);
        assert_eq!(back, text);
        assert!(member.len() <= MOST_BYTES_OF_A_MEMBER, "{}", member.len());
        let size = usize::from(u16::from_le_bytes([member[16], member[17]])) + 1;
        assert_eq!(size, member.len());
    }

    #[test]
    fn a_member_that_stores_its_text_is_one_stored_block_that_decompresses_to_it() {
        let text = bytes_that_do_not_shrink();
        let mut member = HEADER_BEFORE_THE_SIZE.to_vec();
        member.extend_from_slice(&[0, 0]);
        let length = u16::try_from(text.len()).expect("a length");
        store(&[&text], length, &mut member);
        end_the_member(&[&text], length, &mut member).expect("the member");
        // 18 bytes of header, the 5 of a stored block, the text and 8.
        assert_eq!(member.len(), BYTES_OF_THE_HEADER + 5 + TEXT_OF_A_MEMBER + 8);
        assert_eq!(member[BYTES_OF_THE_HEADER], 0x01);
        let mut back = Vec::new();
        MultiGzDecoder::new(&member[..])
            .read_to_end(&mut back)
            .expect("the member decompresses");
        assert_eq!(back, text);
    }

    #[test]
    fn the_deflate_is_kept_up_to_65510_bytes_of_data_of_a_stream_that_ended() {
        assert!(the_deflate_fits(&Ok(Status::StreamEnd), 65510));
        assert!(!the_deflate_fits(&Ok(Status::StreamEnd), 65511));
        assert!(!the_deflate_fits(&Ok(Status::Ok), 100));
        assert!(!the_deflate_fits(&Ok(Status::BufError), 100));
    }

    #[test]
    fn a_member_of_text_that_shrinks_is_deflated_and_states_its_size() {
        let text = b"chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/1\n".repeat(1000);
        let (member, back) = member_and_text(&text[..TEXT_OF_A_MEMBER.min(text.len())]);
        assert_eq!(back, text[..TEXT_OF_A_MEMBER.min(text.len())]);
        assert!(member.len() < 2000, "{}", member.len());
        let size = usize::from(u16::from_le_bytes([member[16], member[17]])) + 1;
        assert_eq!(size, member.len());
    }

    #[test]
    fn a_member_refuses_more_text_than_a_member_holds() {
        let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
        let mut member = Vec::new();
        let text = vec![b'a'; TEXT_OF_A_MEMBER + 1];
        match compress_a_member(&mut compress, &[&text], &mut member) {
            Err(Error::VcfWriterMemberNotBuilt { what, found, most }) => {
                assert_eq!((what, found, most), ("text", 65281, 65280));
            }
            other => panic!("not the error of a member not built: {other:?}"),
        }
    }

    #[test]
    fn the_members_are_cut_across_the_pieces_of_the_text_and_the_tail_waits() {
        let waiting = vec![b'w'; 100];
        let first = vec![b'a'; 65200];
        let second = vec![b'b'; 70000];
        let third = vec![b'c'; 10];
        let (full, tail) = members_of(&waiting, &[&first, &second, &third]);
        let lengths: Vec<Vec<usize>> = full
            .iter()
            .map(|member| member.iter().map(|piece| piece.len()).collect())
            .collect();
        // 100 + 65200 + 70000 + 10 = 135310: two members of 65280 and a
        // tail of 4750.
        assert_eq!(lengths, [vec![100, 65180], vec![20, 65260]]);
        let tail_lengths: Vec<usize> = tail.iter().map(|piece| piece.len()).collect();
        assert_eq!(tail_lengths, [4740, 10]);
        assert_eq!(full[1][0], &first[65180..]);
    }

    #[test]
    fn a_member_of_several_pieces_decompresses_to_their_text_in_order() {
        let pieces: [&[u8]; 3] = [b"chr1\t100\t", b"", b"rs1\tA\tT\n"];
        let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
        let mut member = Vec::new();
        compress_a_member(&mut compress, &pieces, &mut member).expect("the member");
        let mut back = Vec::new();
        MultiGzDecoder::new(&member[..])
            .read_to_end(&mut back)
            .expect("the member decompresses");
        assert_eq!(back, b"chr1\t100\trs1\tA\tT\n");
    }
}
