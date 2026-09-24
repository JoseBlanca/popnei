//! The error of the crate.
//!
//! Every module of popnei adds its cases to [`Error`], and every operation
//! that can fail returns the [`Result`] of this module. The enum is
//! `#[non_exhaustive]`, so a module that is written later adds a case
//! without breaking the code that matches on it.

use std::path::PathBuf;

use thiserror::Error as ThisError;

use crate::block::BlockSize;
use crate::filters::FilteringStats;
use crate::io::vcf::VcfPlace;
use crate::ld::{MAX_ALLELES_OF_A_VARIANT, MAX_PLOIDY_OF_THE_DOSAGES, MAX_VALUES_OF_THE_DOSAGES};
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

    /// The counts of one variant over a population were given the index of
    /// an individual that the variant has no genotype for: the variant
    /// holds one genotype for each individual of the reader, and the index
    /// is at or beyond them.
    ///
    /// The indices of a population are resolved from the names a user wrote
    /// against the individuals the pass gives, before any variant is read,
    /// so a user reaches this only through a defect of popnei.
    #[error(
        "the counts of one variant over a population were given the individual {individual}, and the variant holds the genotypes of {num_individuals} individuals"
    )]
    IndividualBeyondTheVariant {
        /// The index the counts were given.
        individual: usize,
        /// How many individuals the variant holds the genotypes of.
        num_individuals: usize,
    },

    /// A variant that was being turned into one dosage per individual has
    /// more than two different alleles among its called genotypes, and
    /// `transform_to_biallelic` is false. The dosage of a genotype is how
    /// many of its alleles are not the major one, which has a meaning for
    /// two alleles; with the argument true every allele that is not the
    /// major one counts the same. The alleles are those the genotypes
    /// hold and not those the source lists. It is the pass over a row that
    /// refuses it, so every calculation that walks that pass raises it and
    /// takes that argument; the principal components of the variants are
    /// the one that does so far. In Python it is a `ValueError`.
    #[error(
        "the variant at the position {position} among those given has {num_alleles} different alleles among its called genotypes, and the dosage of a genotype, how many of its alleles are not the major one, has a meaning for two: pass `transform_to_biallelic` to count every allele that is not the major one the same"
    )]
    VariantWithMoreThanTwoAlleles {
        /// Which variant of those the reader gave it is, from 0.
        position: usize,
        /// How many different alleles it has among its called genotypes.
        num_alleles: usize,
    },

    /// A variant was to be turned into one dosage per individual at a
    /// ploidy that the dosages cannot be written at. The pass over a row
    /// writes the genotype of each individual as one byte, its dosage or
    /// the code of a genotype with an allele missing, so a ploidy of 255,
    /// the largest the VCF reader takes, has one dosage more than a byte
    /// holds. No organism of this world reaches it: the largest ploidy of
    /// one is a dozen. In Python it is a `ValueError`.
    #[error(
        "a variant cannot be turned into dosages at a ploidy of {ploidy}: the pass writes the genotype of each individual as one byte, its dosage or the code of a genotype with an allele missing, which takes a ploidy of {largest} at most",
        largest = crate::variant::MAX_PLOIDY_OF_THE_VARIANTS
    )]
    VariantPloidyTooLarge {
        /// How many alleles each genotype of the variant holds.
        ploidy: usize,
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

    /// A reader gave a block that holds the genotypes of no individual,
    /// either because it has no individual or because its ploidy is 0, and
    /// a calculation over the variants reads the genotype of one individual
    /// at least. The VCF reader refuses a header with no individual and a
    /// ploidy of 0, so a user reaches this only through a reader with a
    /// defect. It is told apart from the genotypes that nobody asked the
    /// reader for, which are missing from a block for another reason and
    /// leave it empty in the same way.
    #[error(
        "a reader gave a block of {num_individuals} individuals of the ploidy {ploidy}, which holds no genotype for a variant, and a calculation over the variants reads the genotype of 1 individual at least; the reader that gave it has a defect"
    )]
    BlockWithNoGenotypeOfAVariant {
        /// How many individuals the block says it holds the genotypes of.
        num_individuals: usize,
        /// How many alleles the genotype of one individual holds in it.
        ploidy: usize,
    },

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
        /// Which filter it is: `missing_data`, `maf`, `obs_het` or `ld`,
        /// the name its counts have for a Python and a TypeScript user.
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
        /// The kind that is filtered twice: `missing_data`, `maf`,
        /// `obs_het` or `ld`.
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

    /// A filter that compares one number of each variant was asked for with
    /// the criterion of the filter by linkage disequilibrium, which
    /// compares a variant with the variants kept behind it and holds the
    /// window of them between one block and the next. The filter of that
    /// criterion is `LdFilter` and the reader of it is `LdFilteredReader`,
    /// which `chain_of` builds for the criterion, so no call from Python or
    /// from TypeScript reaches this: it is a caller of the core crate that
    /// built the wrong filter for a criterion. In Python it is a
    /// `RuntimeError`, the defect of popnei that it is.
    #[error(
        "a filter that compares one number of each variant was asked for with the criterion of the filter by linkage disequilibrium, which compares a variant with the variants kept behind it: the filter of that criterion is LdFilter, and chain_of builds an LdFilteredReader for it"
    )]
    VarFilterOfTheLdCriterion,

    /// The `max_dist` of the filter by linkage disequilibrium is below 1: a
    /// window that holds nothing but the variants at the very position of
    /// the variant it is the window of. It is the number a user writes, in
    /// `filter_by_ld(0.1, 0)`, so it is refused at the call that adds the
    /// filter, and in Python it is a `ValueError` that names no file. A
    /// negative number never reaches the core: the binding crate takes
    /// `max_dist` as a signed integer and refuses it itself, so that a
    /// user gets that `ValueError` and not the `OverflowError` that pyo3
    /// raises when a negative number is asked of a `u64`.
    #[error(
        "the max_dist of the filter by linkage disequilibrium is {max_dist}, and a window of 0 base pairs reaches no variant but the ones at the very position of the variant it is the window of: a max_dist is 1 base pair or more"
    )]
    LdFilterMaxDistTooSmall {
        /// The distance that was given for it.
        max_dist: u64,
    },

    /// A variant given to the filter by linkage disequilibrium that does
    /// not come after the one before it: its position falls below the
    /// position of the variant before it on its chromosome, or it is on a
    /// chromosome that had already ended. The window of a variant is the
    /// variants kept behind it on its chromosome, so this filter is the one
    /// reader of popnei that needs the variants of each chromosome to come
    /// together and in the order of their positions, where the rest of
    /// popnei reads a source in any order. The alternative is to subtract
    /// two positions that run backwards, which on a `u64` wraps to a
    /// distance of 18 million million million and puts the pair outside
    /// every window without a word. In Python it is a `ValueError` with the
    /// file the variants were read from: it is what that source holds.
    #[error(
        "the variant {variant} of the ones the filter by linkage disequilibrium has read, on {chrom}, does not come after the one before it, and that filter compares a variant with the ones it kept behind it on its chromosome: {problem}; give it a source whose variants come with each chromosome together and in the order of their positions, which `bcftools sort` writes"
    )]
    LdFilterVariantOutOfOrder {
        /// Which variant it is, counted from 1 over the variants the filter
        /// has been given since it was built, which are the ones the
        /// filters before it kept.
        variant: u64,
        /// The chromosome it is on: its name when the reader that gave the
        /// block was there to be asked for it, and the number it has in
        /// that reader's table when the filter was called on a block of its
        /// own.
        chrom: crate::filters::TheChromOfTheVariant,
        /// How it does not come after the variant before it.
        problem: crate::filters::TheOrderOfTheVariants,
    },

    /// A name given to the filter of individuals of
    /// `docs/specs/filters.md` is not an individual of the variants it is
    /// put on. It is the name a user wrote, so the message names it. pyNei
    /// drops it in silence, and `filter_samples(v, ["ind05", "nope"])`
    /// gives variants of one individual.
    #[error(
        "`{name}` is not an individual of the variants; `individuals` gives the names the source has, written as the source writes them"
    )]
    IndividualNotInTheSource {
        /// The name that is not an individual of the variants.
        name: String,
    },

    /// A name given to the filter of individuals is there twice. Two
    /// columns of the genotypes of one individual are one individual for
    /// everything that reads them, and every count over them would hold it
    /// twice. pyNei keeps the individual once.
    #[error(
        "the individual `{name}` is named twice among the individuals to keep, and each of them is kept once"
    )]
    IndividualNamedTwice {
        /// The name that is there twice.
        name: String,
    },

    /// The filter of individuals was given no name at all, which would
    /// leave variants of nobody: every source of popnei holds one
    /// individual at least, as `docs/specs/block.md` says.
    #[error(
        "no individual was named to keep, and the variants hold the genotypes of one individual at least"
    )]
    NoIndividualNamed,

    /// A second filter of individuals on variants that hold one. Two lists
    /// of individuals keep the ones that are in both, which is one list, so
    /// the second says that the user has lost track of the individuals
    /// their variants carry, which running the cell of a notebook twice
    /// gives. pyNei takes it. The case of a second threshold filter is the
    /// one above, which carries the two thresholds that this one has no
    /// counterpart of.
    #[error(
        "the variants are filtered by {kind} already, and a second filter of individuals keeps the individuals that are in both lists, which is one list"
    )]
    FilterOfIndividualsThatIsSet {
        /// The kind of the step, which is `individuals`: the name a Python
        /// and a TypeScript user reads for it.
        kind: &'static str,
    },

    /// A name in one of the populations of `pops` is not an individual of
    /// the variants the statistic is calculated over, which are those of
    /// the source after the filter of individuals when there is one. It is
    /// the name a user wrote, so the message names it and the population it
    /// is in. pyNei refuses it too, in `_calc_pops_idxs`, naming the
    /// population and every name of it that is missing.
    #[error(
        "`{name}` is named in the population `{pop}` and is not an individual of the variants; `individuals` gives the names the variants have, which are the ones the filter of individuals keeps when there is one"
    )]
    IndividualOfAPopNotInThePass {
        /// The population the name was given in.
        pop: String,
        /// The name that is not an individual of the variants.
        name: String,
    },

    /// A name is twice in one population. Every count over the population
    /// would hold that individual twice: pyNei counts it twice, and
    /// `{"x": ["a", "b", "b"]}` over the genotypes `0/0 0/1 1/1` gives it
    /// an observed heterozygosity of 2/3 at commit ef0ca6e. An individual
    /// that is in two populations is taken, as in pyNei.
    #[error(
        "the individual `{name}` is named twice in the population `{pop}`, and a population holds each of its individuals once"
    )]
    IndividualNamedTwiceInAPop {
        /// The population the name is twice in.
        pop: String,
        /// The name that is there twice.
        name: String,
    },

    /// A population of `pops` names no individual. Every statistic of a
    /// population is calculated over its individuals, so a population with
    /// none has no value for any of them; pyNei gives NaN for each.
    #[error(
        "the population `{pop}` names no individual, and every statistic of a population is calculated over the individuals of that population"
    )]
    PopWithNoIndividual {
        /// The population that names no individual.
        pop: String,
    },

    /// `pops` holds no population at all, which would leave a result with
    /// nothing in it: pyNei gives one with no column. A user who wants one
    /// population of every individual gives no `pops`.
    #[error(
        "`pops` names no population, and a result holds one value for each population: leave `pops` out for one population of every individual"
    )]
    NoPop,

    /// The histogram of a statistic was asked for no bin. A histogram
    /// counts the variants that fall in each of its bins, so one with no bin
    /// counts nothing: pyNei hands `num_bins` to `numpy.linspace`, which
    /// gives one edge for 0 bins and a histogram that no value falls in.
    #[error(
        "the histogram was asked for 0 bins, and a histogram has 1 bin at least: `num_bins` is how many bins the values of the statistic are counted in"
    )]
    HistWithNoBin,

    /// The range of the histogram of a statistic does not run from a number
    /// up to a larger one: its two ends are equal, they are the wrong way
    /// round, or one of them is NaN or infinite, which leaves every edge
    /// between them NaN.
    #[error(
        "the range of the histogram is {start:?} to {end:?}, and a range runs from a number up to a larger one: the bins divide that range, and the statistics of one variant lie between 0 and 1"
    )]
    HistRangeNotGoingUp {
        /// The start of the range that was given.
        start: f64,
        /// The end of the range that was given.
        end: f64,
    },

    /// The two ends of the range of a histogram are each a number, and the
    /// distance between them is above the largest float64. The width of a
    /// bin is that distance over the bins, so it is infinite, and the edges
    /// of 4 bins from -1e308 to 1e308 are NaN, infinite, infinite, infinite
    /// and 1e308: they do not go up, and the search for the bin of a value
    /// over edges that do not go up puts every value in the first bin.
    /// numpy refuses the same range, with "Too many bins for data range",
    /// and pyNei's `numpy.histogram` with "'bins' must increase
    /// monotonically".
    #[error(
        "the range of the histogram is {start:?} to {end:?}, and the distance between its two ends is above the largest float64: the width of a bin is that distance over the bins, and the edges of the bins have to go up"
    )]
    HistRangeTooWide {
        /// The start of the range that was given.
        start: f64,
        /// The end of the range that was given.
        end: f64,
    },

    /// The histogram of a statistic was asked for more bins than
    /// `stats::MAX_NUM_BINS`. Every bin is a count of 8 bytes for each
    /// population and each statistic, once in the pass and once more in
    /// every chunk of rows a thread is reading, and a histogram a person
    /// reads has tens of bins.
    #[error(
        "the histogram was asked for {num_bins} bins, and it has {largest} at most: a histogram a person reads has tens of bins, and the counts of more than {largest} of them for each population and each statistic are more memory than a machine gives"
    )]
    HistTooManyBins {
        /// How many bins the histogram was asked for.
        num_bins: usize,
        /// The most it has, `stats::MAX_NUM_BINS`.
        largest: usize,
    },

    /// The range of a histogram whose bins are of equal ratio starts at 0 or
    /// below. Each edge is the one before it times a fixed factor, and no
    /// factor takes 0 anywhere. pyNei refuses it too, in `_prepare_bins`.
    #[error(
        "the range of the histogram starts at {start:?} and its bins are of equal ratio, which start above 0: each edge is the one before it times a fixed factor, and no factor takes 0 anywhere"
    )]
    HistLogRangeNotAboveZero {
        /// The start of the range that was given.
        start: f64,
    },

    /// The ploidy or the exponent a statistic of one variant was built with
    /// is 0, or above the largest ploidy a reader of popnei gives. A trial
    /// implementation of the expected heterozygosity in September 2026 gave,
    /// with an exponent of 0, a plain -2.0 and an unbiased NaN for the
    /// allele counts 2, 1 and 1; a ploidy of 0 turns the
    /// `min_num_individuals` test off, since it asks for 0 called alleles;
    /// and an exponent of 1e8 took 0.2 s for one variant.
    #[error(
        "the {kind} of a statistic of one variant is {value}, and it is 1 at least and {largest} at most, the largest ploidy a reader of popnei gives; the exponent is the number the allele frequencies are raised to, which is the ploidy of the variants unless the caller asks for another one"
    )]
    StatPloidyOutOfRange {
        /// Which of the two it is, `ploidy` or `exponent`.
        kind: &'static str,
        /// The number that was given for it.
        value: usize,
        /// The largest one, `io::vcf::MAX_PLOIDY`.
        largest: usize,
    },

    /// A user asked for a statistic of a variant under a name that is of
    /// none of the five. The names are those of the fields of the result,
    /// and they are `stats::PerVarStat::NAMES`, which the message lists.
    #[error(
        "`{name}` is not one of the statistics of a variant, which are {the_five}",
        the_five = the_five_statistics()
    )]
    StatOfAnUnknownName {
        /// The name the user wrote.
        name: String,
    },

    /// A user asked for bins of a kind that is neither of the two: bins of
    /// equal width, `stats::LINEAR_BINS`, and bins of equal ratio,
    /// `stats::LOGARITHMIC_BINS`. pyNei spells the first one `lineal`, the
    /// Spanish word, and popnei refuses that name as any other unknown one,
    /// which the owner decided on 22 September 2026.
    #[error(
        "`bin_type` is `{kind}`, and the bins of a histogram are `{linear}`, of equal width, or `{logarithmic}`, of equal ratio; pyNei spells the first one `lineal`, the Spanish word",
        linear = crate::stats::LINEAR_BINS,
        logarithmic = crate::stats::LOGARITHMIC_BINS
    )]
    HistBinsOfAnUnknownKind {
        /// The name the user wrote.
        kind: String,
    },

    /// The threshold below which a variant counts as polymorphic in a
    /// population is not a number from 0 to 1. A major allele frequency is
    /// a count of one allele divided by the called alleles, so every value
    /// it takes lies between 0 and 1, and a threshold outside that range
    /// makes every variant polymorphic or none. pyNei compares with
    /// whatever it is given.
    #[error(
        "`poly_threshold` is {value:?}, and it is a number from 0 to 1, both included: a variant is polymorphic in a population when its major allele frequency there is below the threshold, and a frequency lies between 0 and 1"
    )]
    PolyThresholdOutOfRange {
        /// The number that was given for it.
        value: f64,
    },

    /// A pass that calculates a statistic gave no variant, either because
    /// its source holds none or because the steps of the pass kept none of
    /// the variants they were given. A mean over no variant, a histogram
    /// that counts nothing, a rate of an individual over no variant and a
    /// distance between two individuals over no variant say
    /// nothing about a dataset, and a user who gets them has to know which
    /// of the two happened, so the message says it with the variants each
    /// filter was given and kept.
    ///
    /// It is the one case of a pass that gave no variant: every
    /// calculation over a pass raises it, the two of `stats` and
    /// `calc_kosman_sums` of `dists`, and each reads the counts from the
    /// reader it was lent, so that neither binding crate writes the
    /// message again in its own language.
    #[error(
        "{said}",
        said = a_pass_that_gave_no_variant(*num_vars_of_the_source, filters)
    )]
    PassGaveNoVariant {
        /// How many variants the source of the pass gave: what the filter
        /// nearest the source was given, and 0 when the pass has no filter,
        /// since the pass then gave what the source gave.
        num_vars_of_the_source: u64,
        /// Each filter of the chain of the pass with the variants it was
        /// given and the ones it kept, the outermost first, as
        /// `filtering_stats` of a reader gives them.
        filters: Vec<(&'static str, FilteringStats)>,
    },

    /// A value of the table of a principal component analysis is not
    /// finite, an infinity or a NaN, with the place where it is. There is
    /// nothing to give for such a table: the mean of that trait, and with
    /// it every projection, would be a NaN. pyNei refuses a NaN and lets
    /// an infinity through to numpy's decomposition, which raises
    /// `LinAlgError`. In Python it is a `ValueError`.
    #[error(
        "the value at row {row}, trait {col} of the table is {value}, and a principal component analysis needs every value finite"
    )]
    PcaValueNotFinite {
        /// Which row of the table holds it, from 0.
        row: usize,
        /// Which trait of the table holds it, from 0.
        col: usize,
        /// The value that is not finite.
        value: f64,
    },

    /// The table of a principal component analysis is to be standardized
    /// and not centered. Standardizing divides each trait by the standard
    /// deviation it has once it is centered, so the two go together, which
    /// is what pyNei's `do_pca` says as well. In Python it is a
    /// `ValueError`.
    #[error(
        "the table is to be standardized and not centered, and standardizing divides each trait by the standard deviation it has once it is centered: center the table or do not standardize it"
    )]
    PcaStandardizeWithoutCentering,

    /// The table of a principal component analysis has fewer than 2 rows
    /// or no traits. One row has no variation for the components to hold:
    /// pyNei raises the error of the traits with no variance for it when
    /// it standardizes, and without standardizing divides by n - 1 = 0 and
    /// gives percentages that are NaN. In Python it is a `ValueError`.
    #[error(
        "the table is {num_rows} x {num_cols}, and a principal component analysis needs 2 rows at least and 1 trait at least"
    )]
    PcaTableTooSmall {
        /// The rows the table was said to have.
        num_rows: usize,
        /// The traits the table was said to have.
        num_cols: usize,
    },

    /// The traits of a table that is to be standardized and that have no
    /// variance, every value of each one being equal to the others. There
    /// is nothing to divide them by, so the user takes them out or does
    /// not standardize; without standardizing they are no error and get a
    /// weight of 0. In Python it is a `ValueError` whose message names the
    /// traits, since the layer that has the frame puts the name of each
    /// column in the place of its position.
    #[error(
        "{count} of the {num_cols} traits have no variance and cannot be standardized: take them out of the table or do not standardize; they are the traits at {shown}, counting from 0",
        count = positions.len(),
        shown = crate::pca::the_positions_listed(positions)
    )]
    PcaTraitsWithNoVariance {
        /// The position of each trait with no variance among the traits of
        /// the table, from 0 and in order.
        positions: Vec<usize>,
        /// How many traits the table has.
        num_cols: usize,
    },

    /// A trait whose mean or whose standard deviation is not a number the
    /// principal component analysis can use, because the values of that
    /// trait are too large or too small for the arithmetic of an `f64`.
    /// [`crate::pca::TraitScale`] says which of the three it is, and each
    /// of them would otherwise give a result with no meaning: NaN
    /// projections, a trait that quietly becomes a column of zeros, or a
    /// division by 0. The user scales that trait or takes it out. In
    /// Python it is a `ValueError` whose message names the trait, as the
    /// error of the traits with no variance does.
    #[error(
        "the trait at the position {position} cannot be centered or standardized: {problem}; scale that trait or take it out of the table"
    )]
    PcaTraitOutOfRange {
        /// The position of the trait among the traits of the table, from
        /// 0.
        position: usize,
        /// Which of the three it is.
        problem: crate::pca::TraitScale,
    },

    /// No trait of the table of a principal component analysis has
    /// variance once it is centered: every value of every trait is equal
    /// to the others, or the table is all zeros. There is no direction to
    /// give. pyNei gives 0 for every projection and a percentage of NaN
    /// for every component. In Python it is a `ValueError`.
    #[error("no trait has variance, there is nothing to do a PCA with")]
    PcaNoTraitWithVariance,

    /// The buffer of the table of a principal component analysis does not
    /// hold exactly its rows times its traits. Only a caller of the function
    /// of the core crate reaches it, since each binding crate takes the
    /// two numbers from the array it was given, so in Python it is a
    /// `RuntimeError`.
    #[error(
        "the table was given as {num_rows} x {num_cols} and its buffer holds {num_values} values"
    )]
    PcaTableOfAnotherSize {
        /// How many values the buffer holds.
        num_values: usize,
        /// The rows the table was said to have.
        num_rows: usize,
        /// The traits the table was said to have.
        num_cols: usize,
    },

    /// An operation of the crate `popnei-linalg` that a principal
    /// component analysis asked for did not run, with what was being
    /// computed. The dimensions and the values that crate refuses are
    /// checked before it is called, so what is left is a table whose
    /// products are not finite and a machine with too little memory for
    /// the workspace of the eigendecomposition. In Python it is a
    /// `RuntimeError`.
    #[error("the {operation} of the principal component analysis could not be done: {source}")]
    PcaLinalg {
        /// What was being computed: the product of the table with itself,
        /// the product of a block of variants with itself, the
        /// eigendecomposition, or one of the three products that give the
        /// projections of a table, the weights of a table and the weights
        /// of a block of variants.
        operation: &'static str,
        /// What the linear algebra said.
        source: popnei_linalg::Error,
    },

    /// The reader of a principal component analysis of the variants gave
    /// no variant. There is nothing to place the individuals by. It is
    /// pyNei's "There are no variants in the 012 matrix", and in Python it
    /// is a `ValueError`: the steps of the `Variants` let no variant
    /// through, or the source has none.
    #[error("there are no variants to do a PCA with")]
    PcaNoVariants,

    /// No variant of a principal component analysis of the variants has
    /// variance: every one of them has one dosage among its called
    /// genotypes, or no called genotype at all. There is no direction to
    /// give. One individual gives it, since every variant of one
    /// individual has one dosage. In Python it is a `ValueError`.
    ///
    /// The message says the dosage and not the genotype, which is what the
    /// rule reads: a variant whose every genotype is missing, and one of
    /// three alleles read as biallelic whose genotypes are `0/1`, `0/2` and
    /// `0/1`, both reach it with genotypes that differ. pyNei's "Every
    /// variant has the same genotype in every sample" is false of the same
    /// two datasets, and it drops a variant by the same rule.
    /// [`Error::KinshipNoVariantWithVariance`] says the same first half.
    #[error(
        "no variant has more than one dosage among its called genotypes, so none of them varies and there is nothing to do a PCA with"
    )]
    PcaNoVariantWithVariance,

    /// The weights of a principal component analysis of the variants were
    /// asked for and no second pass over the variants was made. The weight
    /// of a variant needs the eigenvectors, which are known when the first
    /// pass ends, so a second reader over the same variants gives them.
    /// Only a caller of the function of the core crate reaches it, since
    /// each binding crate opens both readers, so in Python it is a
    /// `RuntimeError`.
    #[error(
        "the weights of {num_prin_comps} components were asked for and no second pass over the variants was made: the weight of a variant needs the eigenvectors, which are known when the first pass ends, so a second reader over the same variants is given whenever `num_prin_comps` is above 0"
    )]
    PcaSecondPassMissing {
        /// How many components the weights were asked for.
        num_prin_comps: usize,
    },

    /// The second pass of a principal component analysis of the variants
    /// read other variants than the first. The weights it works out belong
    /// to the variants of the first pass, which the eigenvectors come
    /// from, so there is nothing to give. It is what a source that changed
    /// between the two passes gives.
    /// [`crate::pca::VariantsOfTheSecondPass`] says what differed. In
    /// Python it is a `RuntimeError`: no argument is wrong, and the core
    /// has no name of a source to give.
    #[error(
        "the second pass over the variants read other variants than the first: {problem}; the source changed between the two passes"
    )]
    PcaSecondPassDiffers {
        /// What the second pass found that the first did not, or the other
        /// way round.
        problem: crate::pca::VariantsOfTheSecondPass,
    },

    /// The source of a principal component analysis of the variants has no
    /// individual. The components are the axes the individuals of a
    /// dataset are placed on, so there is nobody to place, and the
    /// standardizing of a block would read its rows in chunks of no
    /// allele. Every source of popnei has one individual at least, as
    /// `docs/specs/block.md` says, so it is a caller of the function of the
    /// core crate with a reader of its own that reaches it. In Python it is
    /// a `ValueError`.
    #[error(
        "the source has no individual, and the principal components of the variants are the axes the individuals of a dataset are placed on"
    )]
    PcaNoIndividual,

    /// The second pass of a principal component analysis of the variants
    /// worked out the weight of a variant whose column of the weights is
    /// not there. The variants of that pass are the variants of the first,
    /// which it checks as it goes, so each of them has a column: this is a
    /// defect of popnei, and in Python it is a `RuntimeError`.
    #[error(
        "the weights of the variant at the column {column} of the {num_used} that were used have no column to go in, which is a defect of popnei: the second pass over the variants counts them against the variants of the first and each of them has one"
    )]
    PcaWeightOutOfPlace {
        /// The column the weights were to go in.
        column: usize,
        /// How many variants the first pass used, which is how many
        /// columns the weights have.
        num_used: usize,
    },

    /// A dataset the principal components of its variants cannot be taken
    /// on, because one of its sizes is beyond what the analysis counts in.
    /// [`crate::pca::VariantsTooLarge`] says which of the four it is. In
    /// Python it is a `ValueError`.
    #[error("the principal components of the variants cannot be taken on this dataset: {problem}")]
    PcaVariantsTooLarge {
        /// Which of the four sizes it is, with the number the dataset has.
        problem: crate::pca::VariantsTooLarge,
    },

    /// No variant of a kinship has variance among the individuals it was
    /// asked for: every one of them has one dosage among its called
    /// genotypes, or no called genotype at all. A kinship measures a pair
    /// against the average pair of the panel, and a panel whose variants
    /// give every individual the same dosage has no such average. It is
    /// pyNei's "No variant varies among the samples, there is no kinship",
    /// and in Python it is a `ValueError`.
    ///
    /// The message says the dosage and not the genotype, which is what the
    /// rule reads: a variant where every individual is `0/1`, one where
    /// every genotype is missing, and one of three alleles read as
    /// biallelic where the genotypes are `0/1`, `0/2` and `0/1`, all have
    /// one dosage among their called genotypes and different genotypes.
    #[error(
        "no variant has more than one dosage among its called genotypes, so none of them varies among these individuals and there is no kinship to take"
    )]
    KinshipNoVariantWithVariance,

    /// A value of the matrix of a kinship is not finite, an infinity or a
    /// NaN, with the place where it is. There is nothing to give for such a
    /// matrix: every component would be a NaN. The matrix of a pass is
    /// never one of these, so it is a matrix a user built and then wrote
    /// into, since the checks of the one they build are made when they
    /// build it. In Python it is a `ValueError`.
    ///
    /// The whole matrix is read for it and not the lower half alone, which
    /// is what the components take: a value above the diagonal says the
    /// matrix is wrong as surely as one below it.
    #[error(
        "the value at the row {row}, column {col} of the matrix of the kinship is {value}, and the principal components of a kinship need every value finite"
    )]
    KinshipValueNotFinite {
        /// Which row of the matrix holds it, from 0 among the individuals
        /// of the kinship.
        row: usize,
        /// Which column of it holds it, from 0.
        col: usize,
        /// The value that is not finite.
        value: f64,
    },

    /// A kinship of no individual. A kinship is the matrix of every pair of
    /// a set of individuals, so there is no pair to give. In Python it is a
    /// `ValueError`.
    ///
    /// Three callers reach it. A pass asked for none of the individuals of
    /// its reader, which is an `individuals` of no position; its
    /// components, of a matrix with no row, which is what a user who built
    /// a kinship by hand from an empty frame has; and a pass over a source
    /// that has no individual, which no reader of popnei gives, as
    /// `docs/specs/block.md` says, so that one is a caller of the core
    /// crate with a reader of its own.
    #[error("the kinship has no individual, and a kinship is the matrix of every pair of them")]
    KinshipNoIndividual,

    /// Two individuals of a kinship have no variant called in both of them,
    /// so the sum of their pair would be divided by 0. pyNei divides all
    /// the same and leaves the NaN in the matrix; popnei names the two,
    /// which a user can drop from the panel. The positions are among the
    /// individuals the kinship was asked for, in the order it has them. In
    /// Python it is a `ValueError`.
    ///
    /// The two are the same individual when it has no called genotype at
    /// all among the variants that were used, which is an ordinary
    /// sequencing that failed, and the message then names that one
    /// individual instead of telling a user to drop one of the two.
    #[error(
        "{said}",
        said = a_pair_with_no_variant_called(*one, *other, *num_vars_of_one, *num_vars_of_other)
    )]
    KinshipPairWithNoVariantCalled {
        /// Where the first of the two is among the individuals of the
        /// kinship, from 0.
        one: usize,
        /// Where the second of the two is, from 0.
        other: usize,
        /// How many of the variants that were used are called in the
        /// first.
        num_vars_of_one: u64,
        /// How many of them are called in the second.
        num_vars_of_other: u64,
    },

    /// A dataset a kinship cannot be taken on, because one of its sizes is
    /// beyond what the calculation counts in.
    /// [`crate::kinship::KinshipTooLarge`] says which of the two it is. In
    /// Python it is a `ValueError`.
    #[error("the kinship cannot be taken on this dataset: {problem}")]
    KinshipVariantsTooLarge {
        /// Which of the two sizes it is, with the number the dataset has.
        problem: crate::kinship::KinshipTooLarge,
    },

    /// The linear algebra of a kinship failed. The products of a kinship
    /// are the standardized dosages of a block with themselves and the
    /// genotypes that were called with themselves, both individuals x
    /// individuals, and its principal components are the
    /// eigendecomposition of the matrix. In Python it is a `RuntimeError`:
    /// every size was checked before the work was asked for, so what is
    /// left is a defect of popnei, a matrix whose products are not finite
    /// or a machine with too little memory for the workspace of the
    /// eigendecomposition.
    #[error("the {operation} of the kinship could not be done: {source}")]
    KinshipLinalg {
        /// What was being computed: the product of a block of variants
        /// with itself, the product of the genotypes that were called with
        /// themselves, or the eigendecomposition that gives the principal
        /// components.
        operation: &'static str,
        /// What the linear algebra said.
        source: popnei_linalg::Error,
    },

    /// The distances of that many individuals need more memory than the
    /// machine gives: popnei keeps two `u32` for every pair of them, which
    /// is 8 bytes times the pairs, 400 MB for 10000 individuals, and the
    /// machine did not give them. The individuals are those the reader says
    /// its source has, so the memory is asked for once, when the first
    /// block arrives.
    #[error(
        "the distances of {num_individuals} individuals are {num_pairs} pairs, and this machine did not give the memory of the two counts popnei keeps for each pair, 8 bytes a pair; calculate over fewer individuals"
    )]
    DistancesOfTooManyIndividuals {
        /// How many individuals the source has.
        num_individuals: usize,
        /// How many pairs they make, which a `usize` holds: the
        /// individuals whose pairs are more than one counts are refused by
        /// [`Error::MorePairsThanAreCounted`] before the memory is asked
        /// for.
        num_pairs: usize,
    },

    /// The distances of that many individuals are more pairs than this
    /// machine counts, so popnei cannot give each pair a place, whatever
    /// memory there is: it holds the two counts of the pairs in one vector,
    /// in the order of the distance vector, and a place in a vector is a
    /// `usize`. A `usize` is 64 bits natively and 32 in wasm, where 92682
    /// individuals make 4294930221 pairs and 92683 make more than one
    /// counts.
    #[error(
        "the distances of {num_individuals} individuals are more pairs than this machine counts: popnei gives each pair a place among the others, and a place is counted in a usize, which holds {largest} here; calculate over fewer individuals",
        largest = usize::MAX
    )]
    MorePairsThanAreCounted {
        /// How many individuals the source has.
        num_individuals: usize,
    },

    /// The sums the Kosman distances are worked out from do not fit in a
    /// `u32`. popnei keeps, for each pair of individuals, the ploidy times
    /// the sum of d and how many variants both of them were called at, and
    /// the first is at most the ploidy times the second. So it takes more
    /// than 4295 million variants of the ploidy 1, and 2147 million of the
    /// ploidy 2, in one block or over a whole pass.
    #[error(
        "the Kosman distances of {num_vars} variants of the ploidy {ploidy} add up, for a pair of individuals, beyond the {largest} that popnei keeps for a pair: it keeps the ploidy times the sum of the distances of the pair, which is at most the ploidy times the variants; calculate over fewer variants",
        largest = u32::MAX
    )]
    KosmanSumsTooLarge {
        /// The variants whose distances were being added: the variants of
        /// the block when the sets of bits of that block are built, and the
        /// variants the pass has read so far, that block's among them, when
        /// a block is added to the sums of the pass.
        num_vars: u64,
        /// How many alleles the genotype of one individual holds.
        ploidy: usize,
    },

    /// An index that was given for an individual of a population is not an
    /// individual of the dataset: they are counted from 0, so the last one
    /// of a dataset of n individuals is n − 1. In Python it is a
    /// `ValueError`.
    #[error(
        "the individual {individual} was asked for and the dataset has {num_individuals} individuals, which are counted from 0"
    )]
    LdIndividualNotInTheDataset {
        /// The index that was given.
        individual: usize,
        /// How many individuals the dataset has.
        num_individuals: usize,
    },

    /// The dosages of that many variants of that many individuals are more
    /// values than the linear algebra works a product out over, or than
    /// this machine holds the genotypes of one variant of. The r² of a set
    /// of variants is six products over matrices of the variants times the
    /// individuals, and `crates/popnei-linalg` refuses a matrix of more
    /// values than the routines of BLAS and LAPACK count in, so such
    /// dosages are refused where they are built. In Python it is a
    /// `ValueError`.
    #[error(
        "the dosages of {num_vars} variants of {num_individuals} individuals are more than this machine works r² out over: a matrix of them holds at most {largest} values, which is what the linear algebra counts them in; calculate over fewer variants or over fewer individuals",
        largest = MAX_VALUES_OF_THE_DOSAGES
    )]
    LdDosagesTooLarge {
        /// How many variants the dosages are of.
        num_vars: usize,
        /// How many individuals they were asked for.
        num_individuals: usize,
    },

    /// The variants asked of a set of dosages are not variants of it: a
    /// tile of the products, or a window of the filter by linkage
    /// disequilibrium, that runs past the variants there are. No argument
    /// a user writes asks for a range of variants, so it is a defect of
    /// popnei, and in Python it is a `RuntimeError`.
    #[error("the {asked_for} variants from {first} were asked of dosages of {num_vars} variants")]
    LdRowsNotInTheDosages {
        /// The first variant that was asked for.
        first: usize,
        /// How many variants were asked for from it.
        asked_for: usize,
        /// How many variants the dosages hold.
        num_vars: usize,
    },

    /// The genotypes of the block hold more alleles each than a dosage of
    /// the r² is counted in. No reader of popnei gives such a block: the
    /// VCF reader takes 255 alleles in a genotype at most, and the largest
    /// ploidy of an organism is a dozen. In Python it is a `ValueError`.
    #[error(
        "the genotypes of the block hold {ploidy} alleles each, and a dosage, how many alleles of a genotype are not the major allele of its variant, is counted in one byte, which takes a ploidy of {largest} at most",
        largest = MAX_PLOIDY_OF_THE_DOSAGES
    )]
    LdPloidyTooLarge {
        /// How many alleles the genotype of one individual holds.
        ploidy: usize,
    },

    /// The r² of two sets of dosages that were built over a different
    /// number of individuals was asked for. The sums of a pair run over
    /// the individuals both of its variants were called in, which are the
    /// value at the same place of the two sets, so the two hold the same
    /// individuals of the block in the same order.
    /// [`crate::ld::TheIndividualsThatDiffer`] says how they differ. No
    /// argument a user writes chooses the two sets of one call, so it is a
    /// defect of popnei, and in Python it is a `RuntimeError`.
    #[error(
        "the r² of two sets of dosages built over different individuals was asked for, and the sums of a pair run over the individuals both of its variants were called in, which are the individual at the same place of the two sets: {problem}"
    )]
    LdDosagesOfOtherIndividuals {
        /// How the individuals of the two sets differ.
        problem: crate::ld::TheIndividualsThatDiffer,
    },

    /// The individuals of the block times its ploidy are more alleles in
    /// one variant than the r² comes out of exactly. The four products
    /// the formula takes of the six sums, n·Σxx, n·Σxy, Σx·Σy and (Σx)²,
    /// are each at most the individuals times the ploidy squared, and a
    /// product of two whole numbers is exact in an `f64` while it is at
    /// most 2^53, so above this bound an r² loses digits with nothing to
    /// show for it. No dataset of this world reaches it. In Python it is a
    /// `ValueError` that names no file: calculate over fewer individuals.
    #[error(
        "the {num_individuals} individuals of the block at the ploidy {ploidy} are more than the {largest} alleles of one variant the r² is worked out exactly over, which is 47453132 diploid individuals; calculate over fewer individuals",
        largest = MAX_ALLELES_OF_A_VARIANT
    )]
    LdTooManyAllelesInAVariant {
        /// How many individuals the dosages were to be built over.
        num_individuals: usize,
        /// How many alleles the genotype of one individual holds.
        ploidy: usize,
    },

    /// One of the matrices the r² of a set of variants is worked out
    /// through, and that this machine did not give the memory for: one of
    /// the three matrices of the dosages, or one of the six sums of the
    /// pairs of two sets. The memory is asked for
    /// with `try_reserve_exact`, which gives it back as this error where
    /// `vec![0.0; n]` would end the process, and which also refuses a
    /// matrix whose bytes a `usize` does not count, as one of more than
    /// 2^29 values is in WebAssembly, where a `usize` is 32 bits. In
    /// Python it is a `ValueError` that names no file, since what it
    /// refuses is the size of the calculation and not what any file holds:
    /// calculate over fewer variants or over fewer individuals.
    #[error("this machine has not the memory for {what} of the r², {values} values of 8 bytes")]
    LdNoMemory {
        /// Which matrix could not be allocated, as "How it runs" of
        /// `docs/specs/ld.md` names them.
        what: &'static str,
        /// How many values it holds.
        values: usize,
    },

    /// An individual was asked for more than once when a set of dosages
    /// was built over some of the individuals of a block. A population is
    /// a set of individuals, and one given twice would be counted twice in
    /// the individuals of every pair it is in, in the major allele
    /// frequency of every variant and in each of the six sums. In Python
    /// it is a `ValueError`: it is the indices a user gave for the
    /// individuals of a population.
    #[error(
        "the individual {individual} was asked for more than once, and the individuals of a population are given once each"
    )]
    LdIndividualAskedForTwice {
        /// The index that was given more than once.
        individual: usize,
    },

    /// The buffer given for the r² of two sets of variants does not hold
    /// one value for each pair of them. The caller of the core crate holds
    /// that buffer and no argument a user writes is it, so it is a defect
    /// of popnei, and in Python it is a `RuntimeError`.
    #[error(
        "the r² of {num_vars_of_a} variants against {num_vars_of_b} is one value for each pair of them, and the buffer given holds {num_values}"
    )]
    LdR2OfAnotherSize {
        /// How many values the buffer holds.
        num_values: usize,
        /// How many variants the first set of dosages has.
        num_vars_of_a: usize,
        /// How many variants the second set has.
        num_vars_of_b: usize,
    },

    /// One of the products the r² of two sets of variants is worked out
    /// from that the crate `popnei-linalg` did not do, with the sum it
    /// gives. The three matrices hold whole numbers and their sizes are
    /// checked where the dosages are built, so what is left is a result of
    /// more values than the routines of BLAS and LAPACK count in, the r²
    /// of two sets whose variants multiplied together are more than
    /// 2147483647, which the cap of `calc_r2_matrix` refuses before a user
    /// reaches it. So no argument a user writes is wrong here, and in
    /// Python it is a `RuntimeError`, as the linear algebra of the
    /// principal component analysis is.
    #[error("the {operation} of the r² of two sets of variants could not be worked out: {source}")]
    LdLinalg {
        /// Which of the six sums of the formula was being computed, as
        /// "How it runs" of `docs/specs/ld.md` names them.
        operation: &'static str,
        /// What the linear algebra said.
        source: popnei_linalg::Error,
    },

    /// The pass gave more variants than the matrix of every pair was
    /// allowed to take. The matrix holds one r² for each pair, so it grows
    /// with the square of the variants, 200 MB of `f64` at the 5000 of
    /// [`crate::ld::MAX_NUM_VARS_OF_THE_MATRIX`] and 80 GB at 100000, and
    /// `calc_r2_matrix` stops the pass as soon as it passes the number it
    /// was given rather than reading a dataset it cannot hold the matrix
    /// of. A user who has the memory raises `max_num_vars`, and one who
    /// has not puts a filter on the variants. In Python it is a
    /// `ValueError`.
    #[error(
        "the pass gave {num_vars} variants and `max_num_vars` is {max_num_vars}: the matrix of {num_vars} variants holds one r² for each pair of them, {bytes} bytes of 8 each, and the pass was stopped as soon as it passed that number, so its source may hold more variants; raise `max_num_vars` or filter the variants"
    )]
    LdTooManyVars {
        /// How many variants the pass had given when it was stopped, which
        /// is the first count above `max_num_vars`.
        num_vars: usize,
        /// How many variants the calculation was allowed to take.
        max_num_vars: usize,
        /// How many bytes the matrix of those variants holds, 8 for each
        /// pair.
        bytes: u64,
    },

    /// The `max_num_vars` of the matrix of every pair is more variants
    /// than this machine counts the pairs of: the matrix holds one value
    /// for each pair, which is the variants squared, and that number is
    /// counted in a `usize`, 64 bits natively and 32 in WebAssembly, where
    /// 65536 variants already pass it. The number is looked at before the
    /// pass, so that a cap no matrix could be held under is refused at the
    /// call and not after a dataset has been read. In Python it is a
    /// `ValueError` that names no file.
    #[error(
        "a `max_num_vars` of {max_num_vars} is more variants than the matrix of every pair is counted in: it holds one r² for each pair, which is the variants squared, and this machine counts to {largest}",
        largest = usize::MAX
    )]
    LdMaxNumVarsTooLarge {
        /// How many variants the calculation was allowed to take.
        max_num_vars: usize,
    },
    /// An individual a study was asked to test is not one the source has.
    /// The individuals of a study are given by their position among those
    /// the reader gives, from 0, and this one is at or beyond their count.
    /// The Python and the TypeScript layers turn a name of the phenotype
    /// into that position and are where a name the source has not is
    /// refused, so the core is reached by a caller of `calc_gwas` with
    /// positions of its own. In Python it is a `ValueError`.
    #[error(
        "the individual at the position {individual} was asked to be tested and the source has {num_individuals} individuals, whose positions are 0 to {last}",
        last = num_individuals.saturating_sub(1)
    )]
    GwasIndividualNotInTheDataset {
        /// The position that was asked for, from 0.
        individual: usize,
        /// How many individuals the source has.
        num_individuals: usize,
    },

    /// An individual is twice among the ones a study was asked to test. It
    /// would weigh twice in the null model and in every variant, and its
    /// phenotype would be read at two rows. It is pyNei's repeated
    /// individual of the phenotype, and in Python it is a `ValueError`.
    #[error(
        "the individual at the position {individual} is twice among the ones to test, and each of them is tested once"
    )]
    GwasIndividualTestedTwice {
        /// The position that is there twice, from 0.
        individual: usize,
    },

    /// The individuals a study was asked to test are not in the order the
    /// source has them. Their phenotype, their design and their dosages are
    /// three lists that are read together, row by row, so an order that is
    /// not the source's puts one individual's phenotype against another's
    /// genotypes and the study answers about nobody. The Python and the
    /// TypeScript layers build the positions by walking the individuals of
    /// the source and keeping the ones that have a phenotype, whatever
    /// order the phenotype was given in, so they rise. In Python it is a
    /// `ValueError`.
    #[error(
        "the individual at the position {individual} is to be tested after the one at {after}, and the individuals of a study are tested in the order the source has them, their phenotype and their design in that order too"
    )]
    GwasIndividualsOutOfOrder {
        /// The position that comes too late, from 0.
        individual: usize,
        /// The position it was given after, which is above it.
        after: usize,
    },

    /// A study of no more individuals than its design has columns plus one.
    /// The design holds one column for the intercept and one for each
    /// covariate, the variant adds one more, and what is left over is what
    /// the uncertainty of the variant's effect is measured from: one
    /// individual at least, so the individuals are the columns plus two.
    /// In Python it is a `ValueError`.
    #[error(
        "{num_individuals} individuals are tested and the design has {num_coefs} columns, so the variant would leave nothing to measure its uncertainty from; a study of that design needs the columns plus two individuals"
    )]
    GwasTooFewIndividuals {
        /// How many individuals are tested.
        num_individuals: usize,
        /// How many columns the design has, the intercept among them.
        num_coefs: usize,
    },

    /// A phenotype of a study holds a value that is not finite. The
    /// individuals that are tested are those that have a phenotype, so a
    /// NaN is an individual that should not have been tested at all, and an
    /// infinity would carry through the fit into every variant's effect.
    /// The Python and the TypeScript layers leave out the individuals whose
    /// phenotype is NaN. In Python it is a `ValueError`.
    #[error(
        "the phenotype of the tested individual at the position {position} is {value}, and a study is fitted on numbers; leave that individual out"
    )]
    GwasPhenotypeNotFinite {
        /// Where the value is among the tested individuals, from 0.
        position: usize,
        /// The value that is not finite.
        value: f64,
    },

    /// A phenotype of a binomial trait holds a value that is neither 0 nor 1.
    /// Such a trait is the individuals that have a condition against those
    /// that have not, and a logistic model is fitted to nothing else. It is
    /// pyNei's refusal of a phenotype that is not 0 or 1, and in Python it
    /// is a `ValueError`.
    #[error(
        "the phenotype of the tested individual at the position {position} is {value}, and a binomial trait is 0 or 1"
    )]
    GwasPhenotypeNotBinomial {
        /// Where the value is among the tested individuals, from 0.
        position: usize,
        /// The value that is neither 0 nor 1.
        value: f64,
    },

    /// A value of the design of a study is not a finite number. The Python
    /// and the TypeScript layers refuse a covariate that is missing or is
    /// not a number, so what reaches this is a covariate that came out of
    /// a user's own arithmetic as an infinity, and a caller of the core
    /// crate. Left in, it would reach the rank of the design, which
    /// refuses what it is given, and the user would be told of a defect of
    /// popnei where they gave a wrong covariate. In Python it is a
    /// `ValueError`.
    #[error(
        "the value of the column {coef} of the design at the tested individual {individual} is {value}, and a study is fitted on numbers; the column 0 is the intercept and the others are the covariates in the order they were given"
    )]
    GwasDesignValueNotFinite {
        /// Which tested individual's row it is in, from 0.
        individual: usize,
        /// Which column of the design it is in, from 0, where 0 is the
        /// intercept.
        coef: usize,
        /// The value that is not finite.
        value: f64,
    },

    /// A value of the kinship a study was given is not a finite number. A
    /// mixed model eigendecomposes that matrix before it is fitted, and
    /// one value that is not a number makes every eigenvalue and every
    /// eigenvector one, so the study would come back with a NaN for every
    /// variant and nothing to say which cell it started from.
    /// `Kinship.__post_init__` refuses a matrix that holds a value that is
    /// not a number, so what reaches this is an infinity out of the user's
    /// own arithmetic, and a caller of the core crate. In Python it is a
    /// `ValueError`.
    #[error(
        "the entry of the row {individual} and the column {other} of the kinship is {value}, and a mixed model is fitted on numbers; the rows and the columns are the tested individuals, in the order they were given"
    )]
    GwasKinshipValueNotFinite {
        /// Which tested individual's row it is in, from 0.
        individual: usize,
        /// Which tested individual's column it is in, from 0.
        other: usize,
        /// The value that is not finite.
        value: f64,
    },

    /// Every tested individual of a binomial trait has the same phenotype.
    /// A study of such a trait compares the individuals that have the
    /// condition with those that have not, and one of the two groups is
    /// empty, so no variant can tell them apart. In Python it is a
    /// `ValueError`.
    #[error(
        "every tested individual has the phenotype {value}, and a binomial trait is compared between the individuals that have the condition and those that have not"
    )]
    GwasPhenotypeOfOneValue {
        /// The phenotype they all have, 0 or 1.
        value: f64,
    },

    /// Every tested individual of a continuous trait has the same
    /// phenotype. A study looks for the variants that go with how a trait
    /// differs between the individuals, and a trait that does not differ
    /// has nothing for a variant to go with. Neither of the two things a
    /// fit would do instead is an answer: with no kinship the residual sum
    /// of squares is 0 and the `se` of every variant comes back 0, and
    /// with one the genetic variance is fitted at 0 and the inverse of the
    /// covariance it feeds returns infinities. pyNei refuses a trait of
    /// one value only for a binomial trait, and "Which individuals are
    /// tested, and the design" of `docs/specs/gwas.md` records that
    /// difference. In Python it is a `ValueError`.
    #[error(
        "every tested individual has the phenotype {value}, and a study looks for the variants that go with how a trait differs between the individuals; a trait that is the same in all of them has nothing for a variant to go with"
    )]
    GwasContinuousPhenotypeOfOneValue {
        /// The phenotype they all have.
        value: f64,
    },

    /// The covariates of a study explain the whole of its trait, so the
    /// restricted maximum likelihood of a linear mixed model leaves both
    /// variances at 0, the covariance of the trait is the zero matrix and
    /// its inverse is infinities. It is reached by giving the trait as one
    /// of its own covariates and by covariates that together predict it
    /// exactly, which is a design a user built wrong and not a defect of
    /// popnei: what it gave until 24 September 2026 was the linear
    /// algebra's refusal of a matrix that is not finite, naming an operand
    /// of a product and the file the variants came from, neither of which
    /// is at fault. In Python it is a `ValueError`, as the trait of one
    /// value beside it is.
    #[error(
        "the covariates explain the whole of the trait, so the fit leaves no variance at all: a study looks for the variants that go with what the covariates do not explain, and here there is nothing they do not explain; take out the covariate that carries the trait"
    )]
    GwasDesignExplainsTheTrait,

    /// The kinship a study was given holds two different numbers for one
    /// pair of individuals. A kinship is symmetric, and the
    /// eigendecomposition reads the lower triangle alone, so such a matrix
    /// was being read as its lower half mirrored with no word to the
    /// caller. The `Kinship` of both packages refuses one at the same
    /// tolerance, a share of the largest absolute entry, and this is what
    /// catches a frame written into after it was built and a caller of the
    /// core crate. In Python it is a `ValueError`.
    #[error(
        "the kinship holds {value} for the pair of the tested individuals {individual} and {other} and {and_back} for the same pair the other way round, and a kinship is symmetric; the eigendecomposition reads the lower triangle alone, so the matrix would be read as that half mirrored"
    )]
    GwasKinshipNotSymmetric {
        /// The row of the cell, as a place among the tested individuals.
        individual: usize,
        /// The column of the cell, as such a place.
        other: usize,
        /// What the matrix holds at that cell.
        value: f64,
        /// What it holds at the cell of the same pair the other way round.
        and_back: f64,
    },

    /// The covariance of the working trait of the logistic mixed model,
    /// the kinship times the variance of its random effect plus the
    /// reciprocals of the weights on the diagonal, could not be factored
    /// at the row the value names, counting from 0, so it is not a
    /// covariance. A weight is at most 0.25, so the reciprocals put 4 at
    /// least on every diagonal entry, and what takes such a matrix below 0
    /// is a kinship whose own smallest eigenvalue is below 0 times a
    /// variance large enough to reach it. The per pair denominators of
    /// `docs/specs/kinship.md` are what put that eigenvalue there: a pair
    /// of individuals whose genotypes are missing in different variants is
    /// counted over a different set of variants from the next pair. The
    /// user gives a kinship built from variants with fewer genotypes
    /// missing. It is neither a defect of popnei nor a wrong argument but
    /// the matrix the data made, and in Python it is a `ValueError`.
    #[error(
        "the covariance of the working trait of the logistic mixed model, the kinship times the variance of its random effect plus the weights, could not be factored at its row {at}, counting from 0, so the kinship is not a covariance: missing genotypes leave every pair of individuals counted over its own variants, which can give the matrix an eigenvalue below 0; build the kinship from variants with fewer genotypes missing"
    )]
    GwasKinshipNotACovariance {
        /// The row the factorization stopped at, counting from 0.
        at: usize,
    },

    /// The columns of the design of a study are not independent: a
    /// covariate is constant, or it is a combination of the others, such as
    /// a copy of one or the sum of two. The effects of such a design are
    /// not one set of numbers but many, and the fit would answer with
    /// whichever the arithmetic reached. It is found with the rank of
    /// `popnei-linalg`, how many of the design's columns are independent at
    /// numpy's tolerance, so a design popnei refuses is a design pyNei
    /// refuses. The user takes the covariate out. In Python it is a
    /// `ValueError`.
    #[error(
        "the design has {num_coefs} columns, the intercept among them, and only {rank} of them are independent: a covariate is constant, or it is a combination of the others, such as a copy of one; take it out"
    )]
    GwasCovariatesCollinear {
        /// How many columns the design has, the intercept among them.
        num_coefs: usize,
        /// How many of them are independent.
        rank: usize,
    },

    /// The buffers a study was given do not hold the study it was given:
    /// the phenotype does not hold one value per tested individual, the
    /// design does not hold one row of its columns per tested individual,
    /// or the design has no column at all.
    /// [`crate::gwas::GwasInputShape`] says which of the three it is. Each
    /// binding crate builds the three from the same individuals, so in
    /// Python it is a `RuntimeError`.
    #[error("the study cannot be run on what it was given: {problem}")]
    GwasInputOfAnotherSize {
        /// Which of the three it is, with the sizes that do not agree.
        problem: crate::gwas::GwasInputShape,
    },

    /// An operation of the crate `popnei-linalg` that a study asked for did
    /// not run, with what was being computed. Its dimensions and its values
    /// are checked before it is called: a design of no row or of no column
    /// is [`Error::GwasTooFewIndividuals`] and
    /// [`Error::GwasInputOfAnotherSize`], and a value of it that is not
    /// finite is [`Error::GwasDesignValueNotFinite`]. What is left is a
    /// matrix of more values than that crate takes, which is a design of
    /// more than 2147483647 of them, a machine with too little memory for
    /// the workspace, and a decomposition that did not come out. In Python
    /// it is a `RuntimeError`.
    #[error("the {operation} of the association study could not be done: {source}")]
    GwasLinalg {
        /// What was being computed: the rank of the design, the thin QR a
        /// model is fitted with, a solve against it, or a product of a
        /// block of variants with something the null model holds.
        operation: &'static str,
        /// What the linear algebra said.
        source: popnei_linalg::Error,
    },

    /// The score test was asked of a continuous trait with no kinship.
    /// The only test of a linear model is the t test of the effect it
    /// fitted, and a score test of it would be the same test with the
    /// residual variance held at the null, which no program reports. It is
    /// pyNei's refusal of the same pair, and in Python it is a
    /// `ValueError`.
    #[error(
        "a continuous trait with no kinship is a linear model, whose only test is the t test of the effect it fitted; ask for the Wald test or for none"
    )]
    GwasScoreTestOfALinearModel,

    /// The Wald test was asked of a binomial trait with a kinship. Such a
    /// test fits the model again with each variant in it, and the model
    /// here is a logistic mixed one, so it would be one mixed model fit
    /// for every variant of the dataset. It is pyNei's refusal of the same
    /// pair, and in Python it is a `ValueError`.
    #[error(
        "a binomial trait with a kinship is a logistic mixed model, and a Wald test of it would fit one mixed model for every variant; ask for the score test or for none"
    )]
    GwasWaldTestOfALogisticMixedModel,

    /// The variants a study was given are more than this machine counts
    /// them in, which is 4294967295 in WebAssembly, where a `usize` is 32
    /// bits. Every variant gets a row of the result and is named by its
    /// position among those the reader gave, and neither is a number that
    /// can be counted past the end. In Python it is a `ValueError`.
    #[error(
        "the study was given more variants than this machine counts them in, which is {largest}",
        largest = usize::MAX
    )]
    GwasVariantsTooLarge,

    /// A model of a study answered for another number of variants than the
    /// block it was given holds. It answers for the variants that have
    /// variance among the tested individuals, one `beta`, one `se` and one
    /// `p_value` for each of them, and the message names the column that
    /// is not of that size. It is a defect of popnei, so in Python it is a
    /// `RuntimeError`.
    #[error(
        "the model answered {num_values} values of `{column}` for a block of which {num_with_variance} variants have variance among the tested individuals, and it answers for each of those"
    )]
    GwasAnswersOfAnotherSize {
        /// Which of the three columns is not of that size.
        column: &'static str,
        /// How many values it holds.
        num_values: usize,
        /// How many variants of the block have variance.
        num_with_variance: usize,
    },

    /// The null model of a study was still moving when its fit ended, so
    /// the effects it would report are the ones it happened to be at and
    /// not the ones that fit the trait. The message names the model and
    /// how many rounds it ran.
    ///
    /// A logistic fit reaches it when a covariate separates the
    /// individuals that have the condition from the ones that have not:
    /// there is then no finite effect for that covariate to have, and the
    /// fit walks towards an infinite one. The user takes that covariate
    /// out. Two things end such a fit, and both are the same runaway: the
    /// 50 rounds it is given run out, or the chances it fits reach 0 and 1
    /// and the design weighted by them is no longer a matrix that can be
    /// factored, which stops it earlier. pyNei only meets the first,
    /// because it solves each round with an LU factorization, which
    /// answers a matrix that a Cholesky refuses.
    ///
    /// The logistic mixed model reaches it for a third reason, which is
    /// the kinship, so the message it gets is not the one the two fits
    /// without a kinship get: see
    /// [`the_remedies_of_a_fit_that_did_not_settle`].
    ///
    /// The design alone reaches it too, with another remedy, which is why
    /// the message names two causes. A study whose covariates are so
    /// nearly a combination of each other that the factorization refuses
    /// the system, while the rank check that `Design::of_the_study` makes
    /// with numpy's tolerance lets them through, ends the same way with no
    /// separation anywhere in it. "The logistic model" of
    /// `docs/specs/gwas.md` measures that band on 200 individuals with two
    /// covariates: at a correlation of 1 less 5e-13 the fit runs three
    /// rounds and the factorization refuses the system, and only once the
    /// two covariates agree to within about 1e-14 does the rank check
    /// catch them first. The remedy there is to take one of the two
    /// covariates out, not the one that separates the individuals.
    ///
    /// In Python it is a `ValueError`, as "The logistic mixed model" of
    /// `docs/specs/gwas.md` decides for both fits: pyNei raises a
    /// `RuntimeError` there, and under the rule of `docs/specs/variant.md`
    /// a `RuntimeError` is a defect of popnei where a fit that will not
    /// settle is the data.
    #[error(
        "{what}, and its null model did not settle in the {rounds} rounds it was fitted in. {remedies}",
        what = model.what_it_is_of(),
        remedies = the_remedies_of_a_fit_that_did_not_settle(*model)
    )]
    GwasFitDidNotSettle {
        /// Which of the four models was being fitted, which the message
        /// names with the trait and the kinship that chose it.
        model: crate::gwas::GwasModel,
        /// How many rounds the fit ran before it was given up.
        rounds: usize,
    },

    /// The GRAMMAR-Gamma approximation was asked for by a study with no
    /// kinship. It stands in for the denominator of a mixed model's test,
    /// which is a product with the covariance of the random effect the
    /// kinship is, and a study without one has no such denominator to
    /// approximate. The user gives a kinship or asks for no approximation.
    /// It is pyNei's refusal of the same pair, and in Python it is a
    /// `ValueError`.
    #[error(
        "the GRAMMAR-Gamma approximation stands in for the denominator of a mixed model's test, and a study with no kinship has no such denominator; give a kinship or ask for no approximation"
    )]
    GwasGrammarGammaWithoutAKinship,

    /// The GRAMMAR-Gamma approximation was asked for by a study that has a
    /// kinship, which is the pair it is for, and popnei has not written it
    /// yet. It is refused and not ignored: a study that made the exact test
    /// of every variant and reported that it had approximated nothing would
    /// give the user no way to tell that what they asked for did not
    /// happen. Until it is written the user asks for no approximation and
    /// gets the exact test, which is what every number of
    /// `docs/specs/gwas.md` is. In Python it is a `ValueError`.
    #[error(
        "the GRAMMAR-Gamma approximation is being written; ask for no approximation and every variant gets the exact denominator of its test, which is what it stands in for"
    )]
    GwasGrammarGammaNotBuilt,

    /// The trait and the kinship of a study ask for a model popnei cannot
    /// fit, which the message names. All four are written since 24
    /// September 2026, and what keeps this case is the pass: it chooses a
    /// mixed model only for a study that brought a kinship and then asks
    /// for that kinship again to fit it, and the arm where it is not there
    /// is this error. No study reaches it, since the same pair of the trait
    /// and the kinship chose the model, and it is here rather than a panic
    /// because the core does not panic. In Python it is a `ValueError`,
    /// since it is the study the user asked for that popnei cannot run.
    #[error(
        "popnei cannot run this study yet: {what}, which is being written",
        what = model.what_it_is_of()
    )]
    GwasModelNotBuilt {
        /// Which of the four models the study needs, which the message
        /// names with the trait and the kinship that chose it.
        model: crate::gwas::GwasModel,
    },

    /// A user asked for a trait under a name that is of neither of the
    /// two. The names are `crate::gwas::TraitType::NAMES`, which the
    /// message lists, and both binding crates read them from there. In
    /// Python it is a `ValueError`.
    #[error(
        "`trait` is `{continuous}`, a measurement of each individual, or `{binomial}`, 0 for an individual that has not a condition and 1 for one that has, and `{name}` was given",
        continuous = crate::gwas::TraitType::Continuous.name(),
        binomial = crate::gwas::TraitType::Binomial.name()
    )]
    GwasTraitOfAnUnknownName {
        /// The name that was given.
        name: String,
    },

    /// A user asked for a test under a name that is of neither of the two
    /// popnei makes. The names are `crate::gwas::TestType::NAMES`, which
    /// the message lists. A name that is of a test popnei makes and that
    /// the model of the study has not is another error,
    /// [`Error::GwasScoreTestOfALinearModel`] or
    /// [`Error::GwasWaldTestOfALogisticMixedModel`]. In Python it is a
    /// `ValueError`.
    #[error(
        "`test` is `{wald}`, which fits the model again with the variant in it, or `{score}`, which measures at the null model how steeply the fit would improve if the variant's effect were let off 0, and `{name}` was given",
        wald = crate::gwas::TestType::Wald.name(),
        score = crate::gwas::TestType::Score.name()
    )]
    GwasTestOfAnUnknownName {
        /// The name that was given.
        name: String,
    },

    /// A user asked for a measure of how far apart two populations are
    /// under a name that is of none of the seven. The names are those of
    /// the fields of the result, and they are
    /// `pop_dists::PopDistMeasure::NAMES`, which the message lists. In
    /// Python it is a `ValueError`.
    #[error(
        "`{name}` is not one of the measures of how far apart two populations are, which are {the_seven}",
        the_seven = the_seven_measures()
    )]
    PopDistMeasureOfAnUnknownName {
        /// The name the user wrote.
        name: String,
    },

    /// The variants were asked to be cut into resampling groups of 0 base
    /// pairs. A group is a stretch of one chromosome and holds one base
    /// pair at least. A caller who wants each variant in a group of its own
    /// asks for that, and one who wants no standard error asks for no
    /// groups; neither is a length.
    #[error(
        "the variants were asked to be cut into resampling groups of 0 base pairs, and a group is a stretch of one chromosome 1 base pair long at least"
    )]
    JackknifeGroupOfNoBasePairs,

    /// The distances between populations were asked for fewer than two
    /// populations. Every one of the seven measures is of a pair, so one
    /// population makes no pair and there is nothing to give. In Python it
    /// is a `ValueError`.
    #[error(
        "the distances between populations are calculated for each pair of populations, and `pops` names {num_pops}: name two populations at least"
    )]
    PopDistsOfFewerThanTwoPops {
        /// How many populations the caller named, which is 1: `Pops`
        /// refuses a `pops` that names none.
        num_pops: usize,
    },

    /// The variants of a pass fell into fewer resampling groups than a
    /// standard error is built from. Each group is left out in turn and the
    /// measure calculated again, so a handful of groups gives a number that
    /// says more about where the cuts fell than about the populations, and
    /// a user who chose a length too long for their data is told rather
    /// than handed it. In Python it is a `ValueError`.
    #[error(
        "the variants fell into {num_groups} resampling groups, and a standard error is built from {at_least} at least: cut them into shorter groups, or ask for no standard error"
    )]
    TooFewJackknifeGroups {
        /// How many groups the variants of the pass fell into.
        num_groups: usize,
        /// How many the standard errors need,
        /// [`MIN_NUM_JACKKNIFE_GROUPS`](crate::pop_dists::MIN_NUM_JACKKNIFE_GROUPS).
        at_least: usize,
    },

    /// A source whose variants are cut into resampling groups of a length
    /// gave a variant whose position is below the position of the variant
    /// before it on the same chromosome. The groups are stretches of one
    /// chromosome, cut by comparing the position of a variant with the
    /// first position of the group being filled, so a variant that goes
    /// back joins that group instead of starting one and the groups are
    /// not the stretches the user asked for. "The standard errors" of
    /// `docs/specs/dists.md` has what it does to the standard error. In
    /// Python it is a `ValueError`, and it names the file the variants
    /// were read from.
    #[error(
        "the variants are cut into resampling groups by their position, and the variant at {chrom} {pos} comes after the variant at {chrom} {before} of the same chromosome: a group is a stretch of one chromosome, so sort the source by chromosome and position, or ask for no standard error"
    )]
    JackknifeGroupsVariantGoesBack {
        /// The chromosome of both variants, as the table of the reader
        /// names it, and as its number where that table has no name for
        /// it, which only a reader with a defect gives.
        chrom: String,
        /// The position of the variant that goes back.
        pos: u64,
        /// The position of the variant before it.
        before: u64,
    },

    /// A source whose variants are cut into resampling groups of a length
    /// gave a variant of a chromosome that an earlier variant had left.
    /// The variants of a chromosome that comes back are cut into groups of
    /// their own over the stretch the earlier ones were already cut into,
    /// so the groups overlap and are not the stretches the user asked for.
    /// In Python it is a `ValueError`, and it names the file the variants
    /// were read from.
    #[error(
        "the variants are cut into resampling groups by their position, and the variant at {chrom} {pos} is of a chromosome that the variant at {before_chrom} {before} had left: the variants of one chromosome have to come together, so sort the source by chromosome and position, or ask for no standard error"
    )]
    JackknifeGroupsChromComesBack {
        /// The chromosome that comes back, as the table of the reader
        /// names it.
        chrom: String,
        /// The position of the variant that is on it.
        pos: u64,
        /// The chromosome of the variant before it.
        before_chrom: String,
        /// The position of the variant before it.
        before: u64,
    },

    /// The six sums the distances between populations are worked out from
    /// are more than this machine gave room for: popnei keeps them for each
    /// pair of populations and each resampling group, 48 bytes each, and
    /// either the pairs and the groups are more than a `usize` counts or
    /// the machine did not give their memory. The groups appear while the
    /// variants are read, so it is raised where the sums grow. In Python it
    /// is a `ValueError`.
    #[error(
        "the six sums popnei keeps for each pair of {num_pops} populations within each of the {num_groups} resampling groups the variants have fallen into, 48 bytes each, are more than this machine gave room for: calculate over fewer populations, or cut the variants into longer groups"
    )]
    PopDistSumsTooLarge {
        /// How many populations the pairs are of.
        num_pops: usize,
        /// How many groups the variants read so far have fallen into, and 1
        /// when no standard errors were asked for, since the sums are then
        /// one run of the pairs.
        num_groups: usize,
    },

    /// The populations the distances were asked for make more pairs than
    /// this machine counts, which takes about 93000 of them where a
    /// `usize` is 32 bits, as it is in wasm, and 4294967296 where it is 64.
    /// It is found before a variant is read, so the resampling groups are
    /// none yet and have no part in it, which is what tells it from
    /// [`Error::PopDistSumsTooLarge`]. In Python it is a `ValueError`.
    #[error(
        "{num_pops} populations make more pairs than this machine counts, and every measure of how far apart two populations are is of a pair: calculate over fewer populations"
    )]
    PopDistsOfTooManyPops {
        /// How many populations were given.
        num_pops: usize,
    },

    /// The reader of a pass over the variants says its genotypes hold 0
    /// alleles, or more than the largest ploidy a reader of popnei gives.
    /// The allele frequencies of a population are raised to that ploidy,
    /// and a ploidy of 0 would turn the `min_num_individuals` test off as
    /// well, since it asks for 0 called alleles. The VCF reader refuses
    /// both when it is opened, and the vars file reader refuses a file
    /// whose genotypes hold no allele, so what is left here is a vars file
    /// that says its genotypes hold more alleles than popnei reads. In
    /// Python it is a `ValueError`, and it names the file the variants were
    /// read from.
    #[error(
        "the variants were read at a ploidy of {ploidy}, and the distances between populations are calculated over genotypes of 1 allele at least and {largest} at most: the allele frequencies of a population are raised to the ploidy"
    )]
    PopDistsPloidyOutOfRange {
        /// The ploidy the reader of the pass gives.
        ploidy: usize,
        /// The largest one popnei reads, `io::vcf::MAX_PLOIDY`.
        largest: usize,
    },

    /// The sums of one resampling group do not hold one place for each pair
    /// of the populations the variants are being counted over. Every
    /// variant of a block is added into them pair by pair, so the pairs
    /// after the last place would be counted at no variant while the others
    /// were counted at every one, and each measure of them would come out
    /// of sums of different variants. In Python it is a `RuntimeError`:
    /// nothing a user asks for gives it.
    #[error(
        "the sums of one resampling group hold {num_pairs} pairs, and {num_pops} populations make more; a variant is added into the sums of every pair it counts for"
    )]
    PopDistSumsOfAnotherSize {
        /// How many populations the variants are counted over.
        num_pops: usize,
        /// How many pairs the sums of one group hold.
        num_pairs: usize,
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

/// What [`Error::PassGaveNoVariant`] says: whether the source of the pass
/// held no variant or its steps kept none of the ones they were given, and
/// in the second case what each filter was given and kept.
///
/// `num_vars_of_the_source` is the variants the source gave, and `filters`
/// the counts of the filters of the pass, the outermost first. The counts
/// are said the other way round, the filter nearest the source first, which
/// is the order the variants went through them in.
fn a_pass_that_gave_no_variant(
    num_vars_of_the_source: u64,
    filters: &[(&'static str, FilteringStats)],
) -> String {
    if num_vars_of_the_source == 0 {
        return "the pass gave no variant and its source holds none: a statistic of a \
                pass is calculated over the variants it gives"
            .to_owned();
    }
    let of_each_filter: Vec<String> = filters
        .iter()
        .rev()
        .map(|(kind, stats)| {
            format!(
                "the `{kind}` filter was given {given} and kept {kept}",
                given = stats.vars_processed,
                kept = stats.vars_kept,
            )
        })
        .collect();
    format!(
        "the pass gave no variant: its source gave {num_vars_of_the_source} and the steps \
         kept none of them, {counts}; a statistic of a pass is calculated over the \
         variants it gives",
        counts = of_each_filter.join(", "),
    )
}

/// The seven measures of how far apart two populations are under the names
/// a user writes them, for the message that refuses a name that is of none
/// of them: "`fst`, `f2`, `chord`, `da`, `dest`, `gst` and
/// `gst_standardized`".
fn the_seven_measures() -> String {
    listed(&crate::pop_dists::PopDistMeasure::NAMES)
}

/// What [`Error::KinshipPairWithNoVariantCalled`] says: the two individuals
/// that have no variant called in both of them, or the one individual that
/// has no called genotype at all among the variants that were used.
///
/// The entry of a pair is divided by how many variants both of its
/// individuals were called at, and both cases are that number being 0. A
/// pair reaches it when each of the two was called somewhere and never
/// together; one individual reaches it, against itself, when its sequencing
/// failed, and then every pair it is in has no variant either, so what a
/// user has to do is leave that one out and not one of a pair.
fn a_pair_with_no_variant_called(
    one: usize,
    other: usize,
    num_vars_of_one: u64,
    num_vars_of_other: u64,
) -> String {
    if one == other {
        return format!(
            "the individual at the position {one} has no called genotype among the \
             variants that were used, so its entry of the kinship would be divided by \
             no variant at all; leave it out"
        );
    }
    format!(
        "the individuals at the positions {one} and {other} have no variant called in \
         both of them, so their entry of the kinship would be divided by no variant at \
         all: {num_vars_of_one} {said} called in the first and {num_vars_of_other} in \
         the second; leave one of the two out",
        said = if num_vars_of_one == 1 {
            "variant is"
        } else {
            "variants are"
        },
    )
}

/// The causes of [`Error::GwasFitDidNotSettle`] and the remedy of each,
/// which are the model's own: a fit with a kinship has a third cause that a
/// fit without one has not, and its remedy is the kinship and not a
/// covariate.
///
/// The two covariates of the first message were the whole of it while the
/// logistic regression was the only fit that ran in rounds. The logistic
/// mixed model then took the same message, and a study of one covariate and
/// a kinship that relates every pair alike was told to take one of its two
/// covariates out, which it cannot do and which is not what went wrong.
///
/// A linear model and a linear mixed model are fitted without rounds and
/// reach neither message; they take the one the logistic regression has, so
/// that a case nobody has written yet says something true of any fit.
fn the_remedies_of_a_fit_that_did_not_settle(model: crate::gwas::GwasModel) -> &'static str {
    match model {
        crate::gwas::GwasModel::Lm | crate::gwas::GwasModel::Lmm | crate::gwas::GwasModel::Glm => {
            "Two things do that and they have different remedies: a covariate that separates the individuals that have the condition from the ones that have not has no finite effect for a fit to reach, and the fit walks towards an infinite one, so take that covariate out; or two covariates carry so nearly the same thing that the system of a round can no longer be factored, although they are independent enough for the study to have been accepted, so take one of the two out"
        }
        crate::gwas::GwasModel::Glmm => {
            "Three things do that and they have different remedies: a covariate that separates the individuals that have the condition from the ones that have not has no finite effect for a fit to reach, and the fit walks towards an infinite one, so take that covariate out; or two covariates carry so nearly the same thing that the system of a round can no longer be factored, although they are independent enough for the study to have been accepted, so take one of the two out; or the kinship asks for a random effect that the trait cannot fit, which a kinship that relates every pair alike does, its effect being one number for every individual that the intercept already holds, and then it is the kinship to look at and not a covariate"
        }
    }
}

/// The five statistics of a variant under the names a user writes them, for
/// the message that refuses a name that is of none of them: "`obs_het`,
/// `maf`, `exp_het`, `unbiased_exp_het` and `poly_vars_ratio`".
fn the_five_statistics() -> String {
    listed(&crate::stats::PerVarStat::NAMES)
}

/// `names` in one sentence, each in backticks, the last one after an "and":
/// "`maf` and `obs_het`".
fn listed(names: &[&'static str]) -> String {
    let named: Vec<String> = names.iter().map(|name| format!("`{name}`")).collect();
    match named.split_last() {
        Some((last, before)) => format!("{} and {last}", before.join(", ")),
        // Every table of names this is called with holds names.
        None => String::new(),
    }
}

/// What every operation of popnei that can fail returns.
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::Error;
    use crate::gwas::GwasModel;
    use crate::io::vcf::VcfPlace;
    use crate::variant::Needs;

    /// A fit that did not settle names the causes of the model it was
    /// fitting: the logistic mixed model has the kinship among them, and
    /// the logistic regression, which has no kinship, has not.
    ///
    /// What the message is for is what the user does next, and a study of
    /// one covariate and a kinship that relates every pair alike was being
    /// told to take one of its two covariates out.
    #[test]
    fn the_message_of_a_mixed_fit_that_did_not_settle_names_the_kinship() {
        let of_the_mixed_model = Error::GwasFitDidNotSettle {
            model: GwasModel::Glmm,
            rounds: 200,
        }
        .to_string();
        assert!(
            of_the_mixed_model.contains("it is the kinship to look at"),
            "the logistic mixed model was refused with {of_the_mixed_model}"
        );
        assert!(
            of_the_mixed_model.contains("separates"),
            "the logistic mixed model was refused with {of_the_mixed_model}"
        );
        let of_the_logistic = Error::GwasFitDidNotSettle {
            model: GwasModel::Glm,
            rounds: 50,
        }
        .to_string();
        assert!(
            !of_the_logistic.contains("it is the kinship to look at"),
            "the logistic regression, which has no kinship, was refused with {of_the_logistic}"
        );
        assert!(
            of_the_logistic.contains("separates"),
            "the logistic regression was refused with {of_the_logistic}"
        );
    }

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
