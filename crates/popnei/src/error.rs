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
use crate::variant::{MAX_ALLELE, MISSING_ALLELE, Needs};

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

    /// The alleles given to the counts of one variant are not genotypes
    /// those counts can read: the ploidy is 0, or the alleles are not a
    /// whole number of genotypes of the ploidy.
    ///
    /// A block of a reader of popnei holds, for each variant, one genotype
    /// of the ploidy of the source for each individual, so a user reaches
    /// this only through a reader with a defect.
    #[error(
        "the counts of one variant were given {num_alleles} alleles of the ploidy {ploidy}, and they read one genotype of the ploidy for each individual: the ploidy is 1 at least, and the alleles are a whole number of genotypes of it"
    )]
    GtsNotWholeGenotypes {
        /// How many alleles the counts were given.
        num_alleles: usize,
        /// The ploidy they were to read them as.
        ploidy: usize,
    },

    /// A variant given to the counts of one variant holds more alleles
    /// than a count of them holds. Every count of the two, the genotypes
    /// that are called and the times one allele was seen, is a `u32`, so a
    /// variant of more than 4295 million alleles would be counted into a
    /// number that wrapped.
    ///
    /// One variant holds the individuals of the source times the ploidy
    /// alleles, and no source has that many, so a user reaches this only
    /// through a reader with a defect.
    #[error(
        "the counts of one variant were given {num_alleles} alleles, and a count of the alleles of one variant holds {largest}",
        largest = u32::MAX
    )]
    MoreAllelesThanACountHolds {
        /// How many alleles the counts were given.
        num_alleles: usize,
    },

    /// A genotype given to the counts of one variant holds an allele below
    /// the missing one, -2 or less. Counted as it is, it would be a called
    /// allele of the counts of the genotypes, and it has no place among the
    /// counts of the alleles, which have one entry for each allele from 0
    /// to [`MAX_ALLELE`].
    ///
    /// No reader of popnei gives such an allele, so a user reaches this
    /// only through a reader with a defect. pyNei refuses it too, in
    /// `_count_alleles_per_var`.
    #[error(
        "the genotypes of a variant hold the allele {allele}, and an allele is {MISSING_ALLELE} when it was not called and 0 to {MAX_ALLELE} when it was"
    )]
    AlleleBelowTheMissingOne {
        /// The allele that was found in the genotypes.
        allele: i8,
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

    /// The individuals given to `Block::retain_individuals` hold an index
    /// at or beyond the individuals of the block. They are indices among
    /// the individuals of the block, and the filter of individuals of
    /// `docs/specs/filters.md` gets them from `resolve_individuals`, which
    /// refuses the name that would give one of these, so a user reaches
    /// this only through a reader with a defect. The block is left as it
    /// was.
    #[error(
        "the individuals to keep hold the index {individual}, of a block of {num_individuals} individuals; an individual of a block is an index below how many it holds"
    )]
    IndividualToKeepNotInTheBlock {
        /// The index that is at or beyond the individuals of the block.
        individual: usize,
        /// How many individuals the block holds.
        num_individuals: usize,
    },

    /// The individuals given to `Block::retain_individuals` hold one index
    /// twice, which would leave two columns of the genotypes of one
    /// individual, telling the consumer that there are two individuals
    /// where there is one. `resolve_individuals` refuses the name behind
    /// it, so a user reaches this only through a reader with a defect. The
    /// block is left as it was.
    #[error(
        "the individuals to keep hold the index {individual} twice, and each individual of a block is kept once: two columns of the genotypes of one individual are one individual"
    )]
    IndividualToKeepTwice {
        /// The index that is there twice.
        individual: usize,
    },

    /// `Block::retain_individuals` was given no individual at all, which
    /// would leave a block of nobody's genotypes. Every source of popnei
    /// holds one individual at least, as `docs/specs/block.md` says, and
    /// `resolve_individuals` refuses a filter of individuals that names
    /// none, so a user reaches this only through a reader with a defect.
    /// The block is left as it was.
    #[error(
        "no individual was given to keep of a block, and a block holds the genotypes of one individual at least"
    )]
    NoIndividualToKeep,

    /// The threshold of a filter of variants is not a number from 0 to 1,
    /// both included: it is NaN, it is below 0 or it is above 1. The number
    /// of the variant that the threshold is compared with is one count of
    /// the variant divided by another, so no other threshold says anything
    /// about which variants a user wants.
    ///
    /// It is the number a user writes, in `filter_by_maf(0.95)` and in the
    /// two other methods, so it is refused at the call that adds the
    /// filter. pyNei takes any number: with a negative one no variant
    /// passes, and with one above 1 every variant that has a number does,
    /// so a 95 written for 0.95 filters nothing and says nothing.
    #[error(
        "the threshold of the {kind} filter is {threshold:?}, and a threshold is a number from 0 to 1, both included: the number of the variant it is compared with is one count of the variant divided by another"
    )]
    VarFilterThresholdOutOfRange {
        /// Which filter it is: `missing_data`, `maf` or `obs_het`, the name
        /// its counts have for a Python and a TypeScript user.
        kind: &'static str,
        /// The threshold that was given for it.
        threshold: f64,
    },

    /// A second filter of a kind that the variants are filtered by already.
    /// Two threshold filters of one kind keep the variants that the
    /// stricter of the two keeps alone, so a second one says that the user
    /// has lost track of the filters their variants carry, which running
    /// the cell of a notebook twice gives. pyNei takes it and adds the
    /// counts of the two together.
    #[error(
        "the variants are filtered by {kind} already{set}, and a second filter of that kind, whose threshold is {threshold:?}, keeps the variants that the stricter of the two keeps alone",
        set = threshold_that_is_set
            .map_or_else(String::new, |set| format!(", with a threshold of {set:?}"))
    )]
    VarFilterOfAKindThatIsSet {
        /// The kind that is filtered twice: `missing_data`, `maf` or
        /// `obs_het`.
        kind: &'static str,
        /// The threshold of the filter that was refused, which is the one
        /// the caller wrote.
        threshold: f64,
        /// The threshold of the filter of that kind that is set already,
        /// where whoever refuses the second one knows it. A chain of
        /// readers says which kinds of filter it holds and not with which
        /// thresholds, so a reader built over one gives `None`, and a
        /// binding crate, which has the steps of the variants with their
        /// thresholds, gives the number. The message leaves it out when it
        /// is `None`.
        threshold_that_is_set: Option<f64>,
    },

    /// A step of a pass that `docs/specs/filters.md` describes, that
    /// [`crate::filters::PassStep`] declares and that popnei does not build
    /// yet: the filter of individuals, whose reader and whose method in each
    /// binding crate are still to be written. No user can put such a step
    /// among the steps of their variants, since no method adds one, so
    /// whoever gets this has found a defect of popnei.
    #[error(
        "the `{kind}` step of a pass is declared and not built yet, and popnei was asked for a pass that holds one"
    )]
    PassStepNotBuilt {
        /// The kind of the step, which is `individuals`.
        kind: &'static str,
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

    /// The source the vars file reader was given is not one: it does not
    /// start as an arrow IPC file, its schema has no `popnei` key, the
    /// value of that key is not json or does not hold one of its four
    /// values, or its footer has no `popnei_batches` key. A vars file of
    /// pyNei is refused here, as a file without the `popnei` key.
    #[error("the source is not a vars file: {problem}")]
    NotAVarsFile {
        /// What of a vars file it lacks.
        problem: String,
    },

    /// The vars file says its format is of a version whose first part is
    /// not the one popnei reads. A file whose second part is another one,
    /// a later one too, is read: that is what lets a later version of the
    /// format add a column without making the files or the readers that
    /// are there useless.
    #[error(
        "the vars file says its format is the version {found}, and popnei reads the files whose version starts with {read}; a file of a later version is read by a later popnei"
    )]
    VarsFormatVersion {
        /// The whole version the file gives, `2.0`.
        found: String,
        /// The first part of the version popnei reads, `1`.
        read: &'static str,
    },

    /// A column of the vars file is of another arrow type than the one that
    /// column holds. A column popnei does not know is ignored, as the rule
    /// of the versions asks; one of the six it knows is read as its type
    /// and as no other.
    #[error(
        "the `{column}` column of the vars file is of the arrow type {found}, and popnei reads it as {expected}"
    )]
    VarsColumnType {
        /// The name of the column.
        column: &'static str,
        /// The arrow type it has in the file.
        found: String,
        /// The arrow type popnei reads it as.
        expected: String,
    },

    /// The `gts` column of the vars file holds another number of alleles
    /// for each variant than the individuals and the ploidy of its `popnei`
    /// key give. That width is what turns the flat buffer of the column
    /// into variants, so a file whose two say different things is refused
    /// and not read one allele beside another.
    #[error(
        "the `gts` column of the vars file holds {found} alleles for each variant, and the {num_individuals} individuals of the ploidy {ploidy} that its `popnei` key names hold {expected}"
    )]
    VarsGtsWidth {
        /// How many alleles the column holds for each variant.
        found: usize,
        /// How many the `popnei` key gives, the individuals times the
        /// ploidy.
        expected: usize,
        /// How many individuals that key names.
        num_individuals: usize,
        /// The ploidy it gives.
        ploidy: usize,
    },

    /// A column of the vars file has no value for one of its variants,
    /// where every variant has one: the chromosome, the position, the
    /// alleles and the genotypes. A null `id` is the empty id and a null
    /// `qual` is a variant with no quality, and neither is an error.
    #[error(
        "the `{column}` column of the vars file has no value for its variant {var}, and every variant has one"
    )]
    VarsNullValue {
        /// The name of the column.
        column: &'static str,
        /// Which variant of the file it is, counted from 1, over the whole
        /// file and not inside its batch.
        var: u64,
    },

    /// A variant of the vars file has a quality that is not a finite
    /// number. NaN is what the column of a block holds for a variant with
    /// no quality, which the file writes as a null, so a NaN that is a
    /// value would be read as a variant that has none; and an infinite
    /// quality is a probability of no variant of 0, which is not what phred
    /// scaling says. The VCF reader refuses both for the same reason.
    #[error(
        "the `qual` column of the vars file holds {found} for its variant {var}, and the quality of a variant is a finite number or no value at all, which the file holds as a null"
    )]
    VarsQualityNotFinite {
        /// The value the column holds there.
        found: f32,
        /// Which variant of the file it is, counted from 1, over the whole
        /// file and not inside its batch.
        var: u64,
    },

    /// A genotype of the vars file holds an allele below
    /// [`crate::variant::MISSING_ALLELE`], -1, which is neither an allele
    /// of the variant nor a missing genotype. The alleles of the file are
    /// signed bytes, so a byte of the `gts` column that was damaged after
    /// the file was written, and a file another program wrote, can say one.
    #[error(
        "the `gts` column of the vars file holds the allele {found} for its variant {var}, and an allele is -1, which is the missing one, or a number of 0 or more; a byte of that column was changed after the file was written, and the file has to be fetched or copied again, or the file was written by a program other than popnei, which has to write the alleles of a genotype as -1 and the numbers of the alleles the variant declares"
    )]
    VarsAlleleBelowMissing {
        /// The first allele of the batch that is below the missing one.
        found: i8,
        /// Which variant of the file that allele is of, counted from 1,
        /// over the whole file and not inside its batch.
        var: u64,
    },

    /// The `popnei_batches` key of the footer of the vars file has one
    /// entry for each batch, and this file has another number of one than
    /// of the other, so no entry can be trusted to be that of its batch.
    #[error(
        "the `popnei_batches` key of the footer of the vars file has {found} entries and the file has {expected} batches; there is one entry for each batch"
    )]
    VarsBatchesDoNotMatch {
        /// How many entries the key has.
        found: usize,
        /// How many batches the file has.
        expected: usize,
    },

    /// A batch of the vars file holds another number of variants than its
    /// entry of the footer gives. The number of variants of a file is read
    /// from those entries without reading a batch, so a batch that does not
    /// hold what its entry says is refused when it is read and the number
    /// the file announces is never a wrong one that goes unnoticed.
    #[error(
        "the batch {batch} of the vars file holds {found} variants and its entry of the `popnei_batches` key of the footer says {expected}"
    )]
    VarsBatchNumVars {
        /// Which batch of the file it is, counted from 1. It is a `u64` as
        /// the variant of [`Error::VarsNullValue`] is: both count over the
        /// whole file, which a machine that counts to 4295 million reads
        /// too.
        batch: u64,
        /// How many variants it holds.
        found: usize,
        /// How many its entry says.
        expected: usize,
    },

    /// The buffers of the vars file are compressed with zstd, which no
    /// build of popnei carries: arrow takes zstd from a crate that wraps
    /// the C library, and popnei builds for WebAssembly with no second
    /// compiler. Arrow decompresses a batch when it is read, so this comes
    /// with the first block and not when the file is opened.
    #[error(
        "the buffers of the vars file are compressed with zstd, and popnei reads the files compressed with lz4 and the files with no compression; the file has to be written again with one of those two"
    )]
    VarsZstd,

    /// The `popnei` key of the vars file names one individual twice. The
    /// names are how a user asks for an individual, and two of one name
    /// would be one individual for the user and two columns of genotypes in
    /// the file.
    #[error(
        "the `popnei` key of the vars file names the individual `{name}` twice, and each individual of a file has its own name"
    )]
    VarsIndividualTwice {
        /// The name that is there twice.
        name: String,
    },

    /// A block given to the vars file writer holds another number of
    /// individuals or another ploidy than the writer was built for. The
    /// `popnei` key of the file is written before the first batch, so every
    /// batch holds the individuals that key names. `write_vars` gives the
    /// writer the individuals of its reader, so a user reaches this only
    /// through a reader with a defect.
    #[error(
        "the writer of the vars file was built for {num_individuals} individuals of the ploidy {ploidy} and was given a block of {found_num_individuals} individuals of the ploidy {found_ploidy}; the `popnei` key of a file names the individuals of every one of its batches"
    )]
    VarsBlockDoesNotFit {
        /// How many individuals the writer was built for.
        num_individuals: usize,
        /// The ploidy it was built for.
        ploidy: usize,
        /// How many individuals the block holds.
        found_num_individuals: usize,
        /// The ploidy of the block.
        found_ploidy: usize,
    },

    /// A block given to the vars file writer holds other columns than the
    /// first block it was given. Every batch of an arrow file shares one
    /// schema, and the first block written is what fixes the columns of the
    /// file.
    #[error(
        "the first block written into the vars file holds {first} and a later one holds {found}, so they differ in {differ}; every batch of an arrow file has the columns of one schema, which the first block written fixes",
        differ = first.difference(*found).union(found.difference(*first))
    )]
    VarsBlockColumns {
        /// The fields of the first block written, which are the columns of
        /// the file.
        first: Needs,
        /// The fields of the block that was given now.
        found: Needs,
    },

    /// A block given to the vars file writer holds a chromosome number that
    /// the table of names given with it has no name for. The file holds the
    /// name of the chromosome as text in every row, so a number with no
    /// name cannot be written. The table is the one of the reader the block
    /// came from and has the names of every block that reader gave, so a
    /// user reaches this only through a reader with a defect.
    #[error(
        "a block given to the writer of the vars file holds the chromosome number {number}, and the table of chromosome names given with it has no name for it"
    )]
    VarsChromNameMissing {
        /// The number that has no name.
        number: u32,
    },

    /// A vars file whose genotypes hold no allele: no individual, or the
    /// ploidy 0. It is the writer asked for such a file, and the reader of
    /// one whose `popnei` key names no individual, whose blocks would be
    /// variants of nobody. Every source of popnei has one individual at
    /// least, as `docs/specs/block.md` says, and the `gts` column is in
    /// every vars file.
    #[error(
        "a vars file of {num_individuals} individuals of the ploidy {ploidy} cannot be: it holds the genotypes of one individual at least, of one allele at least each"
    )]
    VarsFileOfNoGenotypes {
        /// How many individuals the writer was asked for.
        num_individuals: usize,
        /// The ploidy it was asked for.
        ploidy: usize,
    },

    /// A column of texts of a block given to the vars file writer holds
    /// more than one column of a batch takes: more bytes, or, in the
    /// `alleles` column, more alleles. Arrow keeps where each text of a
    /// column ends, and where the alleles of each variant end, in a 32 bit
    /// number, and arrow-rs panics at the one that goes past it, so the
    /// writer counts the bytes of the `chrom`, the `id` and the `alleles`
    /// columns of a block, and the alleles of the last, before it fills one
    /// and writes nothing of that block.
    #[error(
        "the `{column}` column of a block given to the writer of the vars file holds {found} {counted}, and one column of a batch of an arrow file holds {largest}; write the file with a smaller `num_vars_per_block`"
    )]
    VarsTextTooLarge {
        /// Which column of the block it is: `chrom`, `id` or `alleles`.
        column: &'static str,
        /// What was counted: `bytes of text`, or `alleles` for the entries
        /// of the `alleles` column.
        counted: &'static str,
        /// How many of those it holds.
        found: u64,
        /// How many one column of a batch holds.
        largest: u64,
    },

    /// The vars file could not be written: the sink refused the bytes, a
    /// disc that filled up among them, or arrow-rs could not write what it
    /// was given.
    ///
    /// It is not [`Error::Io`], which is a source that could not be read. A
    /// call that writes a vars file reads another file, and which of the
    /// two went wrong is what a user acts on, so the write says that it was
    /// the write.
    #[error("the vars file could not be written: {problem}")]
    VarsFileNotWritten {
        /// What went wrong, as the system or arrow-rs said it.
        problem: String,
        /// The error the file system gave, when the cause is one and not a
        /// defect of what arrow-rs was handed. A binding crate builds the
        /// exception of its language with the number it carries, which is
        /// what makes it the `PermissionError` or the `OSError` of that
        /// number in Python.
        source: Option<std::io::Error>,
    },

    /// The vars file starts as an arrow file and ends before what it says
    /// it holds: a download that stopped, a copy that was cut short. The
    /// variants after the cut are not in it, and a reader that gave the
    /// ones before would give fewer variants than the file was written with
    /// and say nothing.
    #[error(
        "the vars file starts as an arrow file and was cut short, so the variants after the cut are not in it and it has to be fetched or copied again: {problem}"
    )]
    VarsFileCutShort {
        /// What was being read when the bytes ran out, as arrow-rs says it.
        problem: String,
    },

    /// A batch of the vars file could not be decoded or decompressed. Its
    /// bytes are not what the file says they are, so the file was damaged
    /// after it was written.
    #[error(
        "the batch {batch} of the vars file could not be read, so the file is damaged and has to be fetched or copied again: {problem}"
    )]
    VarsBatchNotRead {
        /// Which batch of the file it is, counted from 1. It is a `u64` as
        /// the variant of [`Error::VarsNullValue`] is: both count over the
        /// whole file, which a machine that counts to 4295 million reads
        /// too.
        batch: u64,
        /// What arrow-rs said about it.
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

    /// Whoever reports one of these has the genotypes that were counted,
    /// and nothing else: no file and no line, since the counts are given a
    /// row of a block. So the message names the numbers that say which
    /// reader built it wrong.
    #[test]
    fn the_message_of_genotypes_that_are_not_whole_names_the_alleles_and_the_ploidy() {
        let error = Error::GtsNotWholeGenotypes {
            num_alleles: 7,
            ploidy: 2,
        };
        let message = error.to_string();
        assert!(message.contains("7 alleles"), "{message}");
        assert!(message.contains("ploidy 2"), "{message}");

        let error = Error::GtsNotWholeGenotypes {
            num_alleles: 10,
            ploidy: 0,
        };
        let message = error.to_string();
        assert!(message.contains("10 alleles"), "{message}");
        assert!(message.contains("ploidy 0"), "{message}");
        // A ploidy of 0 is read here as one genotype of no allele for
        // every individual, so the message says what a ploidy is.
        assert!(message.contains("the ploidy is 1 at least"), "{message}");
    }

    /// A variant of more alleles than a count of them holds is its own
    /// case, because the counts read every one of those alleles and the
    /// ploidy says nothing about how many there are. The message names
    /// what was given and what a count holds.
    #[test]
    fn the_message_of_more_alleles_than_a_count_holds_names_them_and_the_largest_count() {
        // The alleles are more than a `u32` holds, so the number of the
        // test is written as the largest `usize`, which on the machines
        // popnei builds natively for is 18446744073709551615: a literal
        // above 4295 million is not a `usize` in wasm, where this compiles
        // too.
        let error = Error::MoreAllelesThanACountHolds {
            num_alleles: usize::MAX,
        };
        let message = error.to_string();
        assert!(
            message.contains(&format!("{alleles} alleles", alleles = usize::MAX)),
            "{message}"
        );
        assert!(message.contains("4294967295"), "{message}");
        // It says nothing of a ploidy, which is the case it was taken out
        // of: these alleles are too many whatever the ploidy is.
        assert!(!message.contains("ploidy"), "{message}");
    }

    /// The error of a second filter of one kind is refused in two places,
    /// by a reader over a chain, which knows the kinds of the filters of
    /// that chain and not their thresholds, and by a binding crate, which
    /// has the steps of the variants with both. The message says which
    /// number is which, and says nothing of the threshold that is set when
    /// whoever refused the filter did not have it.
    #[test]
    fn the_message_of_a_second_filter_of_one_kind_names_the_thresholds_it_was_given() {
        let error = Error::VarFilterOfAKindThatIsSet {
            kind: "maf",
            threshold: 0.95,
            threshold_that_is_set: None,
        };
        let message = error.to_string();
        assert!(message.contains("filtered by maf already,"), "{message}");
        assert!(message.contains("whose threshold is 0.95"), "{message}");

        let error = Error::VarFilterOfAKindThatIsSet {
            kind: "maf",
            threshold: 0.95,
            threshold_that_is_set: Some(0.8),
        };
        let message = error.to_string();
        assert!(
            message.contains("filtered by maf already, with a threshold of 0.8"),
            "{message}"
        );
        assert!(message.contains("whose threshold is 0.95"), "{message}");
    }

    /// The allele is what says where the reader that gave it went wrong,
    /// and the range is what says why it was refused.
    #[test]
    fn the_message_of_an_allele_below_the_missing_one_names_the_allele() {
        let error = Error::AlleleBelowTheMissingOne { allele: -2 };
        let message = error.to_string();
        assert!(message.contains("the allele -2"), "{message}");
        assert!(message.contains("-1"), "{message}");
        assert!(message.contains("127"), "{message}");
    }

    /// A user who gets this one has a file that was damaged after it was
    /// written, or one another program wrote, so the message says which
    /// variant of it to look at and what it holds there.
    #[test]
    fn the_message_of_an_allele_of_a_vars_file_below_the_missing_one_names_it_and_its_variant() {
        let error = Error::VarsAlleleBelowMissing { found: -2, var: 17 };
        let message = error.to_string();
        assert!(message.contains("the allele -2"), "{message}");
        assert!(message.contains("variant 17"), "{message}");
        assert!(message.contains("`gts`"), "{message}");
        // A user whose disc changed that byte reads what to do, as they do
        // for the other damaged vars files, and one whose file came from
        // another program reads what that program has to write.
        assert!(message.contains("fetched or copied again"), "{message}");
        assert!(message.contains("a program other than popnei"), "{message}");
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
