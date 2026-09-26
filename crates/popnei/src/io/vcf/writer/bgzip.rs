//! The members of bgzip that a bgzipped VCF is written in, as "How it runs"
//! of the VCF writer in `docs/specs/io_vcf.md` has them: members of 65280
//! bytes of text, each a gzip stream whose header states its size in the
//! extra field `BC`, compressed on the threads of rayon and written in the
//! order of the text, and the empty member that marks the end of the file.

use std::io::Write;

use flate2::{Compress, Compression, Crc, FlushCompress, Status};

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
/// individuals, on one thread of the owner's M5 Pro on 26 September 2026,
/// `write_vcf` bgzipped took 10.3 s at this level and gave 38.4 MB, 4.9 s
/// and 41.6 MB at 5, and 1.7 s and 73.2 MB at 1; `bgzip -@1` took 7.6 s and
/// gave 37.7 MB.
pub(crate) const COMPRESSION_LEVEL: u32 = 6;

/// The first byte of a deflate block that stores its text as it is and is
/// the last block of its stream.
const THE_LAST_STORED_BLOCK: u8 = 0x01;

/// The text of a bgzipped file on its way to the sink: the text that does
/// not fill a member yet waits here from one block to the next, and the
/// members of the rest are compressed and written.
pub(super) struct BgzipOut<W: Write> {
    sink: W,
    /// The text that has not been compressed yet.
    waiting: Vec<u8>,
    /// One buffer for each member compressed at once, kept from one block
    /// to the next.
    members: Vec<Vec<u8>>,
}

impl<W: Write> BgzipOut<W> {
    pub(super) fn new(sink: W) -> BgzipOut<W> {
        BgzipOut {
            sink,
            waiting: Vec::new(),
            members: Vec::new(),
        }
    }

    /// `bytes` of text after the ones given before. They are compressed at
    /// the next [`BgzipOut::write_the_full_members`].
    pub(super) fn write(&mut self, bytes: &[u8]) {
        self.waiting.extend_from_slice(bytes);
    }

    /// The members of every 65280 bytes of the text that waits, compressed
    /// and written, and the rest left waiting.
    ///
    /// # Errors
    ///
    /// When the sink refuses a member, and when a member could not be put
    /// together, which is a defect.
    pub(super) fn write_the_full_members(&mut self) -> Result<()> {
        let num_members = self
            .waiting
            .len()
            .checked_div(TEXT_OF_A_MEMBER)
            .unwrap_or(0);
        if num_members == 0 {
            return Ok(());
        }
        let full = num_members.saturating_mul(TEXT_OF_A_MEMBER);
        let text = self.waiting.get(..full).unwrap_or_default();
        self.members.resize_with(num_members, Vec::new);
        compress_the_members(text, &mut self.members)?;
        for member in &self.members {
            self.sink.write_all(member).map_err(not_written)?;
        }
        self.waiting.drain(..full);
        Ok(())
    }

    /// The member of the text that is left, when there is any, the empty
    /// member of the end, and the sink back with every byte flushed to it.
    ///
    /// # Errors
    ///
    /// When the sink refuses them.
    pub(super) fn finish(mut self) -> Result<W> {
        self.write_the_full_members()?;
        if !self.waiting.is_empty() {
            let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
            let mut member = Vec::new();
            compress_a_member(&mut compress, &self.waiting, &mut member)?;
            self.sink.write_all(&member).map_err(not_written)?;
        }
        self.sink
            .write_all(&THE_EMPTY_MEMBER)
            .map_err(not_written)?;
        self.sink.flush().map_err(not_written)?;
        Ok(self.sink)
    }
}

/// The member of each run of 65280 bytes of `text`, into `members`, one for
/// each run, in their order. Natively the runs are compressed on the
/// threads of rayon, each thread with a deflate of its own.
///
/// # Errors
///
/// What [`compress_a_member`] refuses.
#[cfg(not(target_family = "wasm"))]
fn compress_the_members(text: &[u8], members: &mut [Vec<u8>]) -> Result<()> {
    use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};
    use rayon::slice::ParallelSlice;

    text.par_chunks(TEXT_OF_A_MEMBER)
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
fn compress_the_members(text: &[u8], members: &mut [Vec<u8>]) -> Result<()> {
    let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
    for (text, member) in text.chunks(TEXT_OF_A_MEMBER).zip(members.iter_mut()) {
        compress_a_member(&mut compress, text, member)?;
    }
    Ok(())
}

/// The member of `text`, 65280 bytes at most, into `member`, over what it
/// held: the header with the size of the member, the text deflated, and
/// its CRC32 and length.
///
/// A text that deflate does not shrink into the data a member holds, and a
/// deflate that fails, which miniz_oxide does not do on a buffer it has
/// room in, give a member that stores the text as it is: it always fits,
/// and it decompresses to the same text.
///
/// # Errors
///
/// When the text is longer than a member holds, which is a defect of the
/// caller.
fn compress_a_member(compress: &mut Compress, text: &[u8], member: &mut Vec<u8>) -> Result<()> {
    let too_long = || Error::BlockArrayOfAnotherSize {
        array: "the text of a member of bgzip",
        found: text.len(),
        expected: TEXT_OF_A_MEMBER,
    };
    let length = u16::try_from(text.len())
        .ok()
        .filter(|_| text.len() <= TEXT_OF_A_MEMBER)
        .ok_or_else(too_long)?;
    member.clear();
    member.reserve(MOST_BYTES_OF_A_MEMBER);
    member.extend_from_slice(&HEADER_BEFORE_THE_SIZE);
    member.extend_from_slice(&[0, 0]);
    compress.reset();
    let deflated = compress.compress_vec(text, member, FlushCompress::Finish);
    let data = member.len().saturating_sub(BYTES_OF_THE_HEADER);
    if !matches!(deflated, Ok(Status::StreamEnd)) || data > MOST_DATA_OF_A_MEMBER {
        member.truncate(BYTES_OF_THE_HEADER);
        store(text, length, member);
    }
    end_the_member(text, length, member)
}

/// The data of a member that stores `text`, of `length` bytes, as it is: one
/// deflate block of no compression, after the header in `member`.
fn store(text: &[u8], length: u16, member: &mut Vec<u8>) {
    member.push(THE_LAST_STORED_BLOCK);
    member.extend_from_slice(&length.to_le_bytes());
    member.extend_from_slice(&(!length).to_le_bytes());
    member.extend_from_slice(text);
}

/// The CRC32 and the length of `text` after the data in `member`, and the
/// size of the member in its header.
///
/// # Errors
///
/// When the member is more bytes than its size states in two bytes, which
/// the callers make impossible.
fn end_the_member(text: &[u8], length: u16, member: &mut Vec<u8>) -> Result<()> {
    let mut crc = Crc::new();
    crc.update(text);
    member.extend_from_slice(&crc.sum().to_le_bytes());
    member.extend_from_slice(&u32::from(length).to_le_bytes());
    // The stored member is 65280 bytes of text and 31 of the rest at most,
    // and a deflated one fits by the test above, so its size less 1 fits in
    // two bytes.
    let size = member
        .len()
        .checked_sub(1)
        .and_then(|size| u16::try_from(size).ok())
        .ok_or(Error::BlockArrayOfAnotherSize {
            array: "a member of bgzip",
            found: member.len(),
            expected: MOST_BYTES_OF_A_MEMBER,
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

    use super::{
        BYTES_OF_THE_HEADER, COMPRESSION_LEVEL, HEADER_BEFORE_THE_SIZE, MOST_BYTES_OF_A_MEMBER,
        TEXT_OF_A_MEMBER, compress_a_member, end_the_member, store,
    };

    /// The member of `text`, and the text it decompresses to.
    fn member_and_text(text: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let mut compress = Compress::new(Compression::new(COMPRESSION_LEVEL), false);
        let mut member = Vec::new();
        compress_a_member(&mut compress, text, &mut member).expect("the member");
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
    fn a_member_of_text_that_deflate_cannot_shrink_fits_in_the_bytes_of_a_member() {
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
        store(&text, length, &mut member);
        end_the_member(&text, length, &mut member).expect("the member");
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
        assert!(compress_a_member(&mut compress, &text, &mut member).is_err());
    }
}
