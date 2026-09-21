//! The reader of the members of a source that bgzip wrote.
//!
//! A file made by bgzip is many gzip members one after another, each with
//! 64 KiB of text at most and an empty one at the end, and the header of
//! every member carries the extra field `BC` with the size of that member
//! in the file. [`BgzfReader`] reads such a source by those sizes, as
//! htslib does: it cuts each member from the source and decompresses it on
//! its own, and checks what came out against the CRC32 and the length of
//! the text that the member ends with. A member that is not one bgzip could
//! have written is [`Error::VcfBgzipCorrupted`], and a source that ends
//! before the empty member of the end is [`Error::VcfBgzipEndMissing`].
//!
//! What it is for is that a decoder that goes from one member to the next
//! on its own takes a corrupted file for a whole one: the case, the review
//! that found it and the owner's decision are in "The cases a reader of the
//! rules would not guess" of `docs/specs/io_vcf.md`.
//!
//! Cutting a member and decompressing it are two steps,
//! [`BgzfReader::cut_the_next_member`] and [`decompress_a_whole_member`],
//! and the second is a function that is given the data of one member, a
//! buffer for its text and a decoder, and nothing of the reader: several
//! members can be decompressed at the same time on several threads. Here
//! the two run one after the other, on the thread that reads.
//!
//! What a plan that decompresses side by side has to add: a buffer for the
//! data of every member that is in flight, where this reader reuses one,
//! and a decoder for every thread, where this reader keeps one; and it has
//! to give the text of the members to its caller in the order of the file,
//! which is not the order in which they come out.

use std::io::BufRead;

use flate2::{Crc, Decompress, FlushDecompress, Status};

use crate::error::{Error, Result};

/// The most text a member holds, which is what the length of the text of a
/// member states at most and what bgzip fills one with.
const MOST_TEXT_OF_A_MEMBER: usize = 65536;

/// The room the text of a member is decompressed into: one byte more than
/// the most a member holds, so that a member of exactly that much text is
/// seen to end and one of more is seen not to.
const ROOM_FOR_THE_TEXT: usize = MOST_TEXT_OF_A_MEMBER + 1;

/// The most bytes a member has in the file, which is what the size in its
/// `BC` states at most: the size is held in two bytes as the size less 1.
const MOST_BYTES_OF_A_MEMBER: usize = 65536;

/// The bytes of the header of a gzip member before its extra field: the two
/// of gzip, the method, the flags, the four of the time stamp, the two of
/// the compressor and the system, and the two of the length of the extra
/// field.
const BYTES_BEFORE_THE_EXTRA_FIELD: usize = 12;

/// The bytes after the data of a member: its CRC32 and the length of its
/// text, four each.
const BYTES_AFTER_THE_DATA: usize = 8;

/// The two bytes every gzip member starts with, the method deflate, and the
/// flags of a member of a file that bgzip wrote, which are the one flag that
/// says that the header carries an extra field.
const GZIP_BYTES: [u8; 2] = [0x1f, 0x8b];
const DEFLATE: u8 = 0x08;
const ONLY_AN_EXTRA_FIELD: u8 = 0x04;

/// The two bytes that name the subfield of the extra field in which bgzip
/// writes the size of the member, and the bytes that subfield holds.
const THE_SIZE_SUBFIELD: [u8; 2] = *b"BC";
const BYTES_OF_THE_SIZE: u16 = 2;

/// The bytes of a subfield of an extra field before its own: the two that
/// name it and the two of its length.
const BYTES_BEFORE_A_SUBFIELD: usize = 4;

/// What cutting the next member from the source gave.
enum Member {
    /// The source has no more members: it ended where one member ends and
    /// the next would begin.
    NoMore,
    /// A whole member, with what its last eight bytes state of its text and
    /// how many bytes it holds in the file.
    Whole {
        states: WhatTheMemberStates,
        size: usize,
    },
    /// The source ended inside the member, so what was read of its data is
    /// what there was of it and it has no CRC32 and no length to be checked
    /// against.
    CutShort,
}

/// A reader over the members of a source that bgzip wrote, which gives the
/// text of the file as the members are read.
///
/// It holds one member at a time, its data as the file has it and the text
/// that came out of it, [`MOST_TEXT_OF_A_MEMBER`] each at most, so the
/// memory it uses does not grow with the file.
pub(crate) struct BgzfReader<R: BufRead> {
    source: R,
    /// The deflate data of the member being read, as the file holds it.
    data: Vec<u8>,
    /// The extra field of the header of that member, which holds its size.
    extra_field: Vec<u8>,
    /// The text that came out of it.
    text: Vec<u8>,
    /// How many bytes of that text were consumed.
    consumed: usize,
    /// The decoder of deflate, kept from one member to the next so that
    /// nothing is allocated for each of them.
    inflater: Decompress,
    /// The byte of the source where the member being read starts, which the
    /// error of a member that is corrupted carries.
    offset: u64,
    /// How many members were read whole, so that the one being read is the
    /// one after them.
    members: u64,
    /// Whether the member read last held no text. The empty member is the
    /// mark of the end of a file that bgzip wrote, and bgzip writes one in
    /// the middle of a file too, so what makes it the end is that the
    /// source has no more bytes after it.
    the_last_member_was_empty: bool,
    /// Whether the source ended inside a member or without the empty member
    /// that marks the end of a bgzipped file.
    cut_short: bool,
    /// Whether there is no member left to read.
    done: bool,
}

impl<R: BufRead> BgzfReader<R> {
    /// The reader over `source`, which is the bytes of the file from its
    /// first member on.
    ///
    /// # Errors
    ///
    /// When the machine does not give the 128 KiB the reader holds one
    /// member in.
    pub(crate) fn new(source: R) -> Result<BgzfReader<R>> {
        let mut data = Vec::new();
        room_for(&mut data, MOST_BYTES_OF_A_MEMBER)?;
        let mut text = Vec::new();
        room_for(&mut text, ROOM_FOR_THE_TEXT)?;
        Ok(BgzfReader {
            source,
            data,
            extra_field: Vec::new(),
            text,
            consumed: 0,
            // `false`: the members of a gzip file hold raw deflate data,
            // without the header of zlib.
            inflater: Decompress::new(false),
            offset: 0,
            members: 0,
            the_last_member_was_empty: false,
            cut_short: false,
            done: false,
        })
    }

    /// The text of the source that has not been consumed, which is what is
    /// left of the member being read, and the members after it read one by
    /// one while it is empty. It is empty at the end of the source.
    ///
    /// # Errors
    ///
    /// When a member is not one that bgzip could have written, when the
    /// source cannot be read, and when the machine does not give the memory
    /// of an extra field.
    pub(crate) fn fill(&mut self) -> Result<&[u8]> {
        while !self.done && self.consumed >= self.text.len() {
            if let Err(error) = self.read_the_next_member() {
                // The member whose CRC32 or length is not the one of its
                // text was decompressed before that was known, so the text
                // of it is here: a reader that went on would hand out the
                // very text that the error says is not the text of the
                // file. After an error there is no more.
                self.done = true;
                self.text.clear();
                self.consumed = 0;
                return Err(error);
            }
        }
        Ok(self.text.get(self.consumed..).unwrap_or_default())
    }

    /// `amount` bytes of what [`BgzfReader::fill`] gave are read and are not
    /// given again.
    pub(crate) fn consume(&mut self, amount: usize) {
        self.consumed = self.consumed.saturating_add(amount).min(self.text.len());
    }

    /// The bytes of the next line of the source, with its end of line,
    /// appended to `line`, and how many they were: 0 at the end of the
    /// source.
    ///
    /// # Errors
    ///
    /// The ones of [`BgzfReader::fill`].
    pub(crate) fn read_line(&mut self, line: &mut Vec<u8>) -> Result<usize> {
        let was = line.len();
        let mut read: usize = 0;
        loop {
            let (taken, the_line_ends) = {
                let text = self.fill()?;
                if text.is_empty() {
                    if self.cut_short {
                        // The source ended in the middle of a line: the
                        // bytes of it that arrived are not a line and are
                        // dropped with the lines of the members that are
                        // missing. A whole file whose last line has no end
                        // of line is another thing, and that line is given.
                        line.truncate(was);
                        return Ok(0);
                    }
                    return Ok(read);
                }
                match memchr::memchr(b'\n', text) {
                    Some(at) => {
                        let end = at.saturating_add(1);
                        line.extend_from_slice(text.get(..end).unwrap_or_default());
                        (end, true)
                    }
                    None => {
                        line.extend_from_slice(text);
                        (text.len(), false)
                    }
                }
            };
            self.consume(taken);
            read = read.saturating_add(taken);
            if the_line_ends {
                return Ok(read);
            }
        }
    }

    /// Whether the source ended inside a member or without the empty member
    /// that marks the end of a file that bgzip wrote, which says that the
    /// file was cut short. It is known only when the source has ended.
    pub(crate) fn was_cut_short(&self) -> bool {
        self.cut_short
    }

    /// The next member cut from the source and decompressed into the text,
    /// which is empty when the source has no more members.
    ///
    /// # Errors
    ///
    /// The ones of [`BgzfReader::fill`].
    fn read_the_next_member(&mut self) -> Result<()> {
        self.text.clear();
        self.consumed = 0;
        match self.cut_the_next_member()? {
            Member::NoMore => {
                self.done = true;
                // The file ends with a member that holds no text, and the
                // one read last is that member only when nothing followed
                // it.
                self.cut_short = !self.the_last_member_was_empty;
            }
            Member::CutShort => {
                self.done = true;
                self.cut_short = true;
                self.decompress_what_there_is_of_a_member();
            }
            Member::Whole { states, size } => {
                let which = self.which_member();
                let BgzfReader {
                    data,
                    text,
                    inflater,
                    ..
                } = self;
                decompress_a_whole_member(data, text, inflater, states, which)?;
                self.the_last_member_was_empty = self.text.is_empty();
                self.members = self.members.saturating_add(1);
                self.offset = self
                    .offset
                    .saturating_add(u64::try_from(size).unwrap_or(u64::MAX));
            }
        }
        Ok(())
    }

    /// The next member of the source cut from it by the size its header
    /// states: its data goes into `data` and its extra field into
    /// `extra_field`, and what comes back says whether the member is whole.
    ///
    /// # Errors
    ///
    /// When the header is not one of a member of a file that bgzip wrote,
    /// when the size it states leaves no room for the member, when the
    /// source cannot be read, and when the machine does not give the memory
    /// of the extra field.
    fn cut_the_next_member(&mut self) -> Result<Member> {
        self.data.clear();
        self.extra_field.clear();

        let mut header = [0u8; BYTES_BEFORE_THE_EXTRA_FIELD];
        let read = take_from(&mut self.source, &mut header)?;
        if read == 0 {
            return Ok(Member::NoMore);
        }
        if read < header.len() {
            return Ok(Member::CutShort);
        }
        self.check_the_header(&header)?;

        let bytes_of_the_extra_field = usize::from(two_bytes_of(&header, 10));
        let read = take_from_into(
            &mut self.source,
            &mut self.extra_field,
            bytes_of_the_extra_field,
        )?;
        if read < bytes_of_the_extra_field {
            return Ok(Member::CutShort);
        }
        let size = self.size_of_the_member()?;

        // What is left of the member after its header: its data and the
        // eight bytes that end it. The size holds the whole member, and the
        // header is the bytes before the extra field and the field itself.
        let after_the_header = size
            .checked_sub(BYTES_BEFORE_THE_EXTRA_FIELD)
            .and_then(|left| left.checked_sub(bytes_of_the_extra_field))
            .filter(|left| *left > BYTES_AFTER_THE_DATA)
            .ok_or_else(|| {
                self.corrupted(format!(
                    "the size it states, {size} bytes, leaves no room for its data after the \
                     {BYTES_BEFORE_THE_EXTRA_FIELD} bytes of its header, the \
                     {bytes_of_the_extra_field} of its extra field and the \
                     {BYTES_AFTER_THE_DATA} of its CRC32 and the length of its text"
                ))
            })?;
        let bytes_of_the_data = after_the_header.saturating_sub(BYTES_AFTER_THE_DATA);
        let read = take_from_into(&mut self.source, &mut self.data, bytes_of_the_data)?;
        if read < bytes_of_the_data {
            return Ok(Member::CutShort);
        }
        let mut end = [0u8; BYTES_AFTER_THE_DATA];
        let read = take_from(&mut self.source, &mut end)?;
        if read < end.len() {
            return Ok(Member::CutShort);
        }
        Ok(Member::Whole {
            states: WhatTheMemberStates {
                crc: four_bytes_of(&end, 0),
                text_len: four_bytes_of(&end, 4),
            },
            size,
        })
    }

    /// The first bytes of the header of a member, which say that it is a
    /// gzip member with an extra field and nothing else in its header.
    ///
    /// bgzip writes the flag of the extra field and no other, and BGZF has
    /// the flags of a member fixed to it, so a member with a name or a
    /// comment in its header is one bgzip did not write. Its name would
    /// also move the data of the member, which is cut by a size and not
    /// followed byte by byte.
    ///
    /// # Errors
    ///
    /// When the header is not that of such a member.
    fn check_the_header(&self, header: &[u8]) -> Result<()> {
        if header.get(..GZIP_BYTES.len()) != Some(GZIP_BYTES.as_slice()) {
            return Err(self.corrupted(format!(
                "it starts with {found} and every member of a gzip file starts with `1f 8b`",
                found = crate::io::vcf::shown(header.get(..GZIP_BYTES.len()).unwrap_or_default()),
            )));
        }
        if header.get(2) != Some(&DEFLATE) {
            return Err(self.corrupted(format!(
                "its data is of the method {method} and the data of a member of a bgzip file is \
                 deflate, the method {DEFLATE}",
                method = header.get(2).copied().unwrap_or_default(),
            )));
        }
        if header.get(3) != Some(&ONLY_AN_EXTRA_FIELD) {
            return Err(self.corrupted(format!(
                "its flags are {flags} and the header of a member of a bgzip file carries an \
                 extra field and nothing else, which is the flags {ONLY_AN_EXTRA_FIELD}",
                flags = header.get(3).copied().unwrap_or_default(),
            )));
        }
        Ok(())
    }

    /// The size of the member being read, which its extra field holds in
    /// the subfield `BC`.
    ///
    /// # Errors
    ///
    /// When the subfields do not end where the extra field does, and when
    /// there is no `BC` of two bytes among them.
    fn size_of_the_member(&self) -> Result<usize> {
        let extra_field = the_extra_field(&self.extra_field);
        if let Some(problem) = extra_field.problem {
            return Err(self.corrupted(problem));
        }
        extra_field.size_of_the_member.ok_or_else(|| {
            self.corrupted(format!(
                "its extra field of {bytes} bytes holds no `BC` of two bytes, which is where a \
                 member of a file that bgzip wrote states its size",
                bytes = self.extra_field.len(),
            ))
        })
    }

    /// The data of a member the source ended inside decompressed as far as
    /// it goes, which is the text of the lines that arrived whole before the
    /// cut.
    ///
    /// Nothing here is an error: the member has no CRC32 and no length of
    /// its text to be checked against, since those are among the bytes that
    /// are missing, and what the reader says of such a source is that it was
    /// cut short. A line that the cut fell inside has no end of line, and
    /// the caller drops it: the VCF reader keeps no line of the batch it was
    /// filling when the source ended.
    fn decompress_what_there_is_of_a_member(&mut self) {
        let BgzfReader {
            data,
            text,
            inflater,
            ..
        } = self;
        inflater.reset(false);
        // More bytes of this member could have come, so the decoder is not
        // told that this is the end of its input; what it gives is every
        // byte it could make of the data that arrived. Its error, if the
        // data that arrived is not deflate at all, is the same cut file.
        let _ = inflater.decompress_vec(data, text, FlushDecompress::None);
    }

    /// Which member is being read, which the error of a member that is
    /// corrupted names.
    fn which_member(&self) -> WhichMember {
        WhichMember {
            number: self.members.saturating_add(1),
            offset: self.offset,
        }
    }

    /// The error of the member being read, which names it and the byte of
    /// the source where it starts.
    fn corrupted(&self, problem: String) -> Error {
        corrupted(self.which_member(), problem)
    }
}

/// Which member of a source is being read: its number, counted from 1, and
/// the byte of the source where it starts. It is what the error of a member
/// that is corrupted names it by, and what a decompression that runs away
/// from the reader carries to say which member it was given.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WhichMember {
    pub(crate) number: u64,
    pub(crate) offset: u64,
}

/// What the last eight bytes of a whole member state of its text: its CRC32
/// and how many bytes it holds.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WhatTheMemberStates {
    pub(crate) crc: u32,
    pub(crate) text_len: u32,
}

/// The error of a member that is corrupted, which names it and the byte of
/// the source where it starts.
fn corrupted(which: WhichMember, problem: String) -> Error {
    Error::VcfBgzipCorrupted {
        member: which.number,
        offset: which.offset,
        problem,
    }
}

/// The `data` of a whole member decompressed into `text`, which it empties
/// first, and checked against what the member states of its text.
///
/// It is given the deflate data of one member, a buffer for its text and a
/// decoder, and touches nothing else: the reader that cut the member is not
/// among its arguments, so a caller that cut several members can decompress
/// them at the same time on several threads. `which` is only for the error.
///
/// # Errors
///
/// When the data does not end where the member does, when the decoder
/// refuses it, and when the text that came out is not the text that the
/// CRC32 and the length describe.
fn decompress_a_whole_member(
    data: &[u8],
    text: &mut Vec<u8>,
    inflater: &mut Decompress,
    states: WhatTheMemberStates,
    which: WhichMember,
) -> Result<()> {
    let stated = usize::try_from(states.text_len).unwrap_or(usize::MAX);
    if stated > MOST_TEXT_OF_A_MEMBER {
        return Err(corrupted(
            which,
            format!(
                "it states that its text is {stated} bytes and a member holds \
                 {MOST_TEXT_OF_A_MEMBER} at most"
            ),
        ));
    }
    text.clear();
    inflater.reset(false);
    let status = inflater.decompress_vec(data, text, FlushDecompress::Finish);
    let read = inflater.total_in();
    match status {
        Err(error) => {
            return Err(corrupted(
                which,
                format!("its data could not be decompressed: {error}"),
            ));
        }
        Ok(Status::StreamEnd) if read == u64::try_from(data.len()).unwrap_or(u64::MAX) => {}
        Ok(_) => {
            return Err(corrupted(
                which,
                format!(
                    "its data of {bytes} bytes is not one deflate stream that ends where the \
                     member does: the decoder read {read} of those bytes and gave {out} bytes \
                     of text",
                    bytes = data.len(),
                    out = text.len(),
                ),
            ));
        }
    }
    if text.len() != stated {
        return Err(corrupted(
            which,
            format!(
                "it states that its text is {stated} bytes and the text that came out of it is \
                 {found}",
                found = text.len(),
            ),
        ));
    }
    let mut of_the_text = Crc::new();
    of_the_text.update(text);
    if of_the_text.sum() != states.crc {
        return Err(corrupted(
            which,
            format!(
                "it states the CRC32 {crc:08x} of its text and the text that came out of it has \
                 the CRC32 {found:08x}",
                crc = states.crc,
                found = of_the_text.sum(),
            ),
        ));
    }
    Ok(())
}

/// What the extra field of the header of a gzip member holds of what makes
/// a member of a file that bgzip wrote.
pub(crate) struct TheExtraField {
    /// The size of the member, which the subfield `BC` of two bytes states
    /// when the field carries one. It is what says that bgzip wrote the
    /// member, since nothing else writes that subfield.
    pub(crate) size_of_the_member: Option<usize>,
    /// Why the subfields are not the subfields of an extra field: one of
    /// them ends after the field does, or the field ends in the middle of
    /// one. A member of a bgzip file whose extra field has this is
    /// corrupted; the first member of a source is asked only whether it has
    /// a `BC`, so that a source whose members are to be checked is read by
    /// the reader that checks them.
    pub(crate) problem: Option<String>,
}

/// The subfields of an extra field walked: BGZF lets a member carry others
/// beside the `BC`, before it or after it, so the field is read subfield by
/// subfield and not at a fixed place. Each of them names itself in two
/// bytes, gives its length in two more and holds that many.
pub(crate) fn the_extra_field(extra_field: &[u8]) -> TheExtraField {
    let of_the_field = extra_field.len();
    let mut at: usize = 0;
    while let Some(left) = of_the_field.checked_sub(at).filter(|left| *left > 0) {
        if left < BYTES_BEFORE_A_SUBFIELD {
            return TheExtraField {
                size_of_the_member: None,
                problem: Some(format!(
                    "the subfields of its extra field of {of_the_field} bytes leave {left} bytes \
                     over at its end, and a subfield names itself in two bytes and gives its \
                     length in two more"
                )),
            };
        }
        let of_the_subfield = two_bytes_of(extra_field, at.saturating_add(2));
        let data = at.saturating_add(BYTES_BEFORE_A_SUBFIELD);
        let end = data.saturating_add(usize::from(of_the_subfield));
        if end > of_the_field {
            return TheExtraField {
                size_of_the_member: None,
                problem: Some(format!(
                    "a subfield of its extra field of {of_the_field} bytes says that it holds \
                     {of_the_subfield} bytes and ends {over} bytes after the field does",
                    over = end.saturating_sub(of_the_field),
                )),
            };
        }
        let names_the_size =
            extra_field.get(at..at.saturating_add(2)) == Some(THE_SIZE_SUBFIELD.as_slice());
        if names_the_size && of_the_subfield == BYTES_OF_THE_SIZE {
            // `BC` holds the size of the whole member less 1, so it fits in
            // two bytes.
            return TheExtraField {
                size_of_the_member: Some(
                    usize::from(two_bytes_of(extra_field, data)).saturating_add(1),
                ),
                problem: None,
            };
        }
        at = end;
    }
    TheExtraField {
        size_of_the_member: None,
        problem: None,
    }
}

/// `out.len()` bytes of `source`, or every byte left when it holds fewer,
/// and how many they were.
///
/// # Errors
///
/// When the source cannot be read.
fn take_from<R: BufRead>(source: &mut R, out: &mut [u8]) -> Result<usize> {
    let mut filled: usize = 0;
    while filled < out.len() {
        let buffer = match source.fill_buf() {
            Ok(buffer) => buffer,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(Error::Io(error)),
        };
        if buffer.is_empty() {
            break;
        }
        let end = out.len().min(filled.saturating_add(buffer.len()));
        let Some(into) = out.get_mut(filled..end) else {
            break;
        };
        let taken = into.len();
        into.copy_from_slice(buffer.get(..taken).unwrap_or_default());
        source.consume(taken);
        filled = end;
    }
    Ok(filled)
}

/// `amount` bytes of `source` appended to `out`, or every byte left when it
/// holds fewer, and how many they were.
///
/// # Errors
///
/// When the source cannot be read and when the machine does not give the
/// memory. `amount` is at most the bytes of one member, which the caller
/// has checked.
fn take_from_into<R: BufRead>(source: &mut R, out: &mut Vec<u8>, amount: usize) -> Result<usize> {
    room_for(out, amount)?;
    let mut filled: usize = 0;
    while filled < amount {
        let buffer = match source.fill_buf() {
            Ok(buffer) => buffer,
            Err(error) if error.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(Error::Io(error)),
        };
        if buffer.is_empty() {
            break;
        }
        let taken = buffer.len().min(amount.saturating_sub(filled));
        out.extend_from_slice(buffer.get(..taken).unwrap_or_default());
        source.consume(taken);
        filled = filled.saturating_add(taken);
    }
    Ok(filled)
}

/// Room for `bytes` in `out`, which holds nothing.
///
/// The memory is asked for with `try_reserve`, which gives it back as an
/// error where `Vec::with_capacity` ends the process. Every size that
/// reaches it is a size the file states, and is at most the 65536 bytes of
/// one member.
///
/// # Errors
///
/// When the machine does not give the memory, as an error of the input:
/// nothing of the file is wrong with it.
fn room_for(out: &mut Vec<u8>, bytes: usize) -> Result<()> {
    out.try_reserve(bytes).map_err(|_| {
        Error::Io(std::io::Error::new(
            std::io::ErrorKind::OutOfMemory,
            format!("the {bytes} bytes of a member of the bgzipped source were not given"),
        ))
    })
}

/// The two bytes of `bytes` at `at`, which every number of a gzip header is
/// held in, the smaller first.
fn two_bytes_of(bytes: &[u8], at: usize) -> u16 {
    let of_the_number = bytes
        .get(at..at.saturating_add(2))
        .and_then(|two| <[u8; 2]>::try_from(two).ok())
        .unwrap_or([0, 0]);
    u16::from_le_bytes(of_the_number)
}

/// The four bytes of `bytes` at `at`, which the CRC32 and the length of the
/// text of a member are held in, the smaller first.
fn four_bytes_of(bytes: &[u8], at: usize) -> u32 {
    let of_the_number = bytes
        .get(at..at.saturating_add(4))
        .and_then(|four| <[u8; 4]>::try_from(four).ok())
        .unwrap_or([0, 0, 0, 0]);
    u32::from_le_bytes(of_the_number)
}

#[cfg(test)]
pub(crate) mod tests {
    use std::io::Cursor;

    use super::{
        BgzfReader, Decompress, ROOM_FOR_THE_TEXT, WhatTheMemberStates, WhichMember,
        decompress_a_whole_member, four_bytes_of, room_for,
    };
    use crate::error::Error;

    /// The empty gzip member of 28 bytes that bgzip writes at the end of a
    /// file, which says that the file is whole: a source that bgzip wrote
    /// and that does not end with a member that holds no text was cut
    /// short. htslib writes these same bytes in every file it closes, and
    /// the tests that write a file that bgzip could have written end it
    /// with them.
    ///
    /// The reader compares nothing with these bytes: what it asks of a file
    /// is that its last member hold no text, which the CRC32 and the length
    /// of that member say.
    pub(crate) const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];

    /// The bytes of a member of a file that bgzip wrote: `text` compressed
    /// with raw deflate, a header with the extra field `BC` that holds the
    /// size of the whole member less 1, and the CRC32 of the text and its
    /// length after the data.
    ///
    /// `before_bc` are the bytes of another subfield of the extra field,
    /// written before `BC`, which BGZF allows and bgzip does not write.
    pub(crate) fn bgzf_member_with(text: &[u8], before_bc: &[u8]) -> Vec<u8> {
        let mut encoder =
            flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::fast());
        std::io::Write::write_all(&mut encoder, text).unwrap();
        let data = encoder.finish().unwrap();

        let extra_field = before_bc.len().checked_add(6).unwrap();
        let total = extra_field
            .checked_add(12)
            .and_then(|bytes| bytes.checked_add(data.len()))
            .and_then(|bytes| bytes.checked_add(8))
            .unwrap();
        assert!(total <= 65536, "a member holds 65536 bytes at most");
        let mut member = vec![0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(&u16::try_from(extra_field).unwrap().to_le_bytes());
        member.extend_from_slice(before_bc);
        member.extend_from_slice(b"BC");
        member.extend_from_slice(&2u16.to_le_bytes());
        member.extend_from_slice(
            &u16::try_from(total.checked_sub(1).unwrap())
                .unwrap()
                .to_le_bytes(),
        );
        member.extend_from_slice(&data);
        let mut crc = flate2::Crc::new();
        crc.update(text);
        member.extend_from_slice(&crc.sum().to_le_bytes());
        member.extend_from_slice(&u32::try_from(text.len()).unwrap().to_le_bytes());
        assert_eq!(member.len(), total);
        member
    }

    /// A member as bgzip writes one, whose extra field is the `BC` alone.
    pub(crate) fn bgzf_member(text: &[u8]) -> Vec<u8> {
        bgzf_member_with(text, &[])
    }

    /// A file that bgzip could have written: one member for the text of
    /// each of `texts` and the empty member of 28 bytes that marks the end.
    pub(crate) fn bgzf_file(texts: &[&[u8]]) -> Vec<u8> {
        let mut file = Vec::new();
        for text in texts {
            file.extend_from_slice(&bgzf_member(text));
        }
        file.extend_from_slice(&BGZF_EOF);
        file
    }

    /// A reader that gave an error gives no more text: the member whose
    /// CRC32 failed is decompressed when the error is found, and a reader
    /// that went on would hand out the text that the error says is not the
    /// text of the file. The VCF reader stops at the first error it gets,
    /// so nothing of popnei asks for a block after one; this is the
    /// contract of the module, which the plan that decompresses several
    /// members at a time builds on.
    #[test]
    fn after_an_error_the_reader_of_the_members_gives_no_text_of_the_members_after_it() {
        let mut file = bgzf_member(b"one\n");
        let mut second = bgzf_member(b"two\n");
        // The CRC32 of the second member, which are the four bytes before
        // the length of its text.
        let at = second.len().checked_sub(8).unwrap();
        second[at] ^= 0x01;
        file.extend_from_slice(&second);
        file.extend_from_slice(&bgzf_member(b"three\n"));
        file.extend_from_slice(&BGZF_EOF);

        let mut reader = BgzfReader::new(Cursor::new(file)).expect("the reader");
        let mut line = Vec::new();
        assert_eq!(reader.read_line(&mut line).expect("the first line"), 4);
        assert_eq!(line, b"one\n");

        let error = reader
            .read_line(&mut line)
            .expect_err("the error of the CRC32");
        let Error::VcfBgzipCorrupted { member, .. } = error else {
            panic!("the error is {error}");
        };
        assert_eq!(member, 2);

        for call in 0..3 {
            line.clear();
            assert_eq!(
                reader.read_line(&mut line).expect("no more"),
                0,
                "call {call}"
            );
            assert!(line.is_empty(), "call {call}: {line:?}");
            assert!(reader.fill().expect("no more").is_empty(), "call {call}");
        }
        // The source was not cut short: its bytes are all there and one of
        // its members is not what it says it is.
        assert!(!reader.was_cut_short());
    }

    /// The decompression of a member does not touch the reader that cut it:
    /// the data of two members, cut as the reader cuts them, are
    /// decompressed on two threads at once, each with its own buffer and
    /// its own decoder. It is what the plan that decompresses the members
    /// of one file side by side will do, and the module doc says what that
    /// plan still has to add.
    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn the_data_of_two_members_is_decompressed_on_two_threads_at_once() {
        let members = [bgzf_member(b"one\n"), bgzf_member(b"two\n")];
        let texts = std::thread::scope(|threads| {
            let working: Vec<_> = members
                .iter()
                .enumerate()
                .map(|(number, member)| {
                    threads.spawn(move || {
                        // The header of a member bgzip writes is 18 bytes
                        // and its last eight are its CRC32 and the length
                        // of its text.
                        let data = member.get(18..member.len() - 8).unwrap();
                        let states = WhatTheMemberStates {
                            crc: four_bytes_of(member, member.len() - 8),
                            text_len: four_bytes_of(member, member.len() - 4),
                        };
                        let which = WhichMember {
                            number: number as u64 + 1,
                            offset: 0,
                        };
                        let mut text = Vec::new();
                        room_for(&mut text, ROOM_FOR_THE_TEXT).unwrap();
                        let mut inflater = Decompress::new(false);
                        decompress_a_whole_member(data, &mut text, &mut inflater, states, which)
                            .expect("the text of the member");
                        text
                    })
                })
                .collect();
            working
                .into_iter()
                .map(|working| working.join().expect("the thread"))
                .collect::<Vec<Vec<u8>>>()
        });
        assert_eq!(texts, [b"one\n".to_vec(), b"two\n".to_vec()]);
    }
}
