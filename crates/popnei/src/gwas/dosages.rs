//! The dosages of a block: for each variant, how many copies of an allele
//! that is not the major one each tested individual carries.
//!
//! [`GwasDosages`] holds them for one block, with the frequency of each
//! variant among the tested individuals and whether it has any variance
//! there, and [`GwasDosages::read_the_block`] is what a pass calls for each
//! block. Every model of the study is fitted and tested on what it gives.

use std::num::NonZeroUsize;

use crate::block::Block;
use crate::error::{Error, Result};
use crate::variant::{
    AlleleCounts, MAX_PLOIDY_OF_THE_VARIANTS, MISSING_CODE, Needs, count_alleles,
    the_codes_of_the_genotypes, the_major_allele, the_row_of,
};

use super::study::{Design, GwasInput};

/// Whether a variant with more than two alleles among the called genotypes
/// of the tested individuals is read or refused, which
/// [`GwasInput::transform_to_biallelic`] chooses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum MultiallelicVariants {
    /// Such a variant is an error naming its position among the variants
    /// the reader has given.
    Refused,
    /// Every allele that is not the major one counts the same, so the
    /// dosage of a genotype is how many of its alleles are not the major
    /// one, whichever they are.
    Collapsed,
}

impl MultiallelicVariants {
    /// What the study asked for.
    pub(super) fn of_the_study(input: &GwasInput<'_>) -> MultiallelicVariants {
        match input.transform_to_biallelic {
            true => MultiallelicVariants::Collapsed,
            false => MultiallelicVariants::Refused,
        }
    }
}

/// What reading one variant gives beside the dosages themselves.
#[derive(Debug, Clone, Copy)]
struct RowDosages {
    /// The frequency of the alleles that are not the major one among the
    /// tested individuals, which is the mean dosage over the ploidy.
    allele_freq: f64,
    /// Whether the called genotypes of the tested individuals hold two
    /// dosages at least, which is what a variant needs to be tested.
    has_variance: bool,
}

/// The buffers one thread keeps while it reads the dosages of the rows of
/// a block, so that nothing is allocated for a variant.
struct DosageScratch {
    /// How often each allele was called among the tested individuals,
    /// which gives the major allele of the variant and how many different
    /// alleles it has.
    allele_counts: AlleleCounts,
    /// The code of the genotype of each tested individual: its dosage, 0
    /// to the ploidy, or [`MISSING_CODE`] for a genotype with an allele
    /// missing.
    codes: Vec<u8>,
}

impl DosageScratch {
    /// The buffers of one thread, for the rows of `num_individuals` tested
    /// individuals.
    fn of(num_individuals: usize) -> DosageScratch {
        DosageScratch {
            allele_counts: [0; 128],
            codes: vec![0; num_individuals],
        }
    }
}

/// The dosages of one variant over the tested individuals written into
/// `row`, with the frequency of the alleles that are not the major one and
/// whether the variant has any variance among them.
///
/// `gts` is the genotypes of one variant, `ploidy` alleles for each tested
/// individual and nobody else, and `row` holds one value for each of them.
/// `position` is which variant of those the reader gave this one is, which
/// the error of a variant with more than two alleles names.
///
/// The dosage of a genotype is how many of its alleles are not the major
/// one, and a genotype with any allele missing takes the mean dosage of
/// its variant, so that once the variant is centered it pulls the
/// individual in no direction. Both the major allele and that mean are of
/// the tested individuals alone, as "What it gives" of
/// `docs/specs/gwas.md` asks, and so is the frequency. The dosages are not
/// divided by anything: `beta` of the result is the effect of one more
/// copy of a non major allele, in the units of the trait, and dividing
/// them by the deviation of the variant would give it in deviations
/// instead. That is where this differs from the standardized row of
/// [`crate::variant`], which every variant of a kinship or of a principal
/// component analysis goes through, and it is why the mean comes back
/// here: the study reports it, as `allele_freq`.
///
/// A variant with no called genotype among the tested individuals has a
/// mean of nothing, which is 0, as `_calc_dosages` of pyNei sets it: its
/// dosages are all 0, its frequency is 0 and it has no variance.
///
/// # Errors
///
/// [`Error::VariantPloidyTooLarge`] when `ploidy` is above
/// [`MAX_PLOIDY_OF_THE_VARIANTS`], which the dosages could not be written
/// one to a byte at. [`Error::GtsNotWholeGenotypes`] when `ploidy` is 0 or
/// `gts` does not hold one genotype of it for each value of `row`.
/// [`Error::VariantWithMoreThanTwoAlleles`] when the tested individuals
/// have more than two different alleles among their called genotypes and
/// the study did not ask for those variants to be read. And whatever the
/// counts of the alleles of one variant refuse, which is
/// [`Error::AlleleBelowTheMissingOne`] and
/// [`Error::MoreAllelesThanACountHolds`].
fn the_dosages_of_a_row(
    gts: &[i8],
    ploidy: usize,
    position: usize,
    multiallelic: MultiallelicVariants,
    scratch: &mut DosageScratch,
    row: &mut [f64],
) -> Result<RowDosages> {
    if ploidy > MAX_PLOIDY_OF_THE_VARIANTS {
        return Err(Error::VariantPloidyTooLarge { ploidy });
    }
    let num_individuals = row.len();
    // One genotype of the ploidy for each tested individual. The ploidy of
    // 0 that `NonZeroUsize` refuses is among these: it would be a genotype
    // of no allele for every one of them.
    let of_a_genotype = match NonZeroUsize::new(ploidy) {
        Some(of_a_genotype) if gts.len() == num_individuals.saturating_mul(ploidy) => of_a_genotype,
        _ => {
            return Err(Error::GtsNotWholeGenotypes {
                num_alleles: gts.len(),
                ploidy,
            });
        }
    };
    let DosageScratch {
        allele_counts,
        codes,
    } = scratch;
    // The alleles are counted for the major allele and for how many
    // different alleles the variant has. How many of them were called is
    // not what the mean is taken over: a genotype with one allele called
    // and one missing has a called allele and no dosage.
    count_alleles(gts, allele_counts)?;
    let num_alleles = allele_counts.iter().filter(|count| **count > 0).count();
    match multiallelic {
        MultiallelicVariants::Refused if num_alleles > 2 => {
            return Err(Error::VariantWithMoreThanTwoAlleles {
                position,
                num_alleles,
            });
        }
        MultiallelicVariants::Refused | MultiallelicVariants::Collapsed => {}
    }
    // The buffer of the codes belongs to the thread and is as long as the
    // rows it has read so far, which are all of the tested individuals:
    // this asks the machine for nothing after the first row.
    codes.resize(num_individuals, 0);
    the_codes_of_the_genotypes(gts, of_a_genotype, the_major_allele(allele_counts), codes);
    let called = the_called_dosages(codes);
    let mean = match called.genotypes {
        0 => 0.0,
        genotypes => called.dosages as f64 / genotypes as f64,
    };
    for (value, code) in row.iter_mut().zip(codes.iter().copied()) {
        *value = match code {
            MISSING_CODE => mean,
            dosage => f64::from(dosage),
        };
    }
    Ok(RowDosages {
        allele_freq: mean / ploidy as f64,
        has_variance: called.highest > called.lowest,
    })
}

/// What the called genotypes of one variant hold: how many they are, the
/// sum of their dosages, and the lowest and the highest of those dosages.
#[derive(Debug, Clone, Copy)]
struct CalledDosages {
    /// How many genotypes were called.
    genotypes: u64,
    /// The sum of their dosages.
    dosages: u64,
    /// The lowest dosage among them, and [`u8::MAX`] when none was called.
    lowest: u8,
    /// The highest dosage among them, and 0 when none was called.
    highest: u8,
}

/// What the called genotypes of one variant hold, read from the code of
/// each of the tested individuals.
///
/// A genotype with an allele missing has no dosage: it is counted in none
/// of the four, so the mean is of the called genotypes and the two
/// extremes are theirs. A variant with no called genotype ends with a
/// lowest of [`u8::MAX`] and a highest of 0, which is why a variant has
/// variance at the strict comparison of the two and not at an inequality.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "a dosage is at most the ploidy, 254, and the genotypes of a variant are at \
              most the alleles of it, which the counts of those alleles checked to be a \
              number a u32 holds, so the sum is below 2^40"
)]
fn the_called_dosages(codes: &[u8]) -> CalledDosages {
    let mut called = CalledDosages {
        genotypes: 0,
        dosages: 0,
        lowest: u8::MAX,
        highest: 0,
    };
    for code in codes.iter().copied() {
        let missing = code == MISSING_CODE;
        called.genotypes += u64::from(!missing);
        called.dosages += u64::from(if missing { 0 } else { code });
        called.lowest = called.lowest.min(if missing { u8::MAX } else { code });
        called.highest = called.highest.max(if missing { 0 } else { code });
    }
    called
}

/// The dosages of the rows of a block over the tested individuals written
/// into `dosages`, and what each row gave, in the order of the block.
///
/// The rows are read on the threads of rayon, as section 3 of
/// `docs/architecture.md` asks: no row reads another and each one writes
/// its own values, so neither the dosages nor the frequencies depend on
/// how many threads there are. The threads are those of the pool the
/// caller is running in, and rayon's global pool only when the caller is
/// in none. Each thread keeps the buffers of one row and allocates nothing
/// per variant.
///
/// `gts` holds the rows of the block, `alleles_per_var` alleles each, over
/// the tested individuals alone; `dosages` holds one value for each of
/// those individuals for each of those rows. `first_var` is which variant
/// of those the reader has given the first row of the block is.
///
/// The error is the one of the first row of the block that has one,
/// wherever the threads found it: each row gives its own result and they
/// are read in the order of the block, so a user who reports a file gets
/// the same message every time.
///
/// # Errors
///
/// What reading the dosages of one row refuses, and
/// [`Error::GwasVariantsTooLarge`] when the position of a variant is
/// beyond what a `usize` counts.
#[cfg(not(target_family = "wasm"))]
fn the_dosages_of_the_rows(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    use rayon::slice::{ParallelSlice, ParallelSliceMut};

    let rows: Vec<Result<RowDosages>> = gts
        .par_chunks_exact(alleles_per_var.get())
        .zip(dosages.par_chunks_exact_mut(num_individuals.get()))
        .enumerate()
        .map_init(
            || DosageScratch::of(num_individuals.get()),
            |scratch, (var, (gts, row))| {
                let position = first_var
                    .checked_add(var)
                    .ok_or(Error::GwasVariantsTooLarge)?;
                the_dosages_of_a_row(gts, ploidy, position, multiallelic, scratch, row)
            },
        )
        .collect();
    rows.into_iter().collect()
}

/// The same rows, read one after another, which is what WebAssembly does:
/// it has no threads.
///
/// # Errors
///
/// The same as the rows read on threads.
#[cfg(target_family = "wasm")]
fn the_dosages_of_the_rows(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    the_dosages_of_the_rows_one_by_one(
        gts,
        alleles_per_var,
        num_individuals,
        ploidy,
        multiallelic,
        first_var,
        dosages,
    )
}

/// The rows read one after another into the buffers of one row: what
/// WebAssembly does, and what the test that compares the two ways of
/// reading a block calls.
///
/// # Errors
///
/// The same as the rows read on threads.
#[cfg(any(target_family = "wasm", test))]
fn the_dosages_of_the_rows_one_by_one(
    gts: &[i8],
    alleles_per_var: NonZeroUsize,
    num_individuals: NonZeroUsize,
    ploidy: usize,
    multiallelic: MultiallelicVariants,
    first_var: usize,
    dosages: &mut [f64],
) -> Result<Vec<RowDosages>> {
    let mut scratch = DosageScratch::of(num_individuals.get());
    let mut rows = Vec::new();
    for (var, (gts, row)) in gts
        .chunks_exact(alleles_per_var.get())
        .zip(dosages.chunks_exact_mut(num_individuals.get()))
        .enumerate()
    {
        let position = first_var
            .checked_add(var)
            .ok_or(Error::GwasVariantsTooLarge)?;
        rows.push(the_dosages_of_a_row(
            gts,
            ploidy,
            position,
            multiallelic,
            &mut scratch,
            row,
        )?);
    }
    Ok(rows)
}

/// What the pass over the blocks knows about the block it is handing over:
/// how many alleles the genotype of one individual holds, which the reader
/// gives and which every block of a dataset has, and which variant of
/// those the reader has given the first row of this block is.
///
/// The two are a struct and not two arguments because they are both a
/// count of something and a caller that swapped them would compile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BlockOfThePass {
    /// How many alleles the genotype of one individual holds, which the
    /// reader says its source has. A block of another ploidy is refused:
    /// the frequency of a variant is its mean dosage over this number, so
    /// a block read at another one gives frequencies of nothing.
    pub ploidy: usize,
    /// Which variant of those the reader has given the first row of the
    /// block is, counted from 0. It is what the error of a variant with
    /// more than two alleles names, so it is the reader's count and not
    /// the block's.
    pub first_var: usize,
}

/// The dosages of the variants of one block over the individuals a study
/// tests, with the frequency of the alleles that are not the major one of
/// each variant and whether it has any variance, both over those
/// individuals alone.
///
/// A dosage is how many alleles of a genotype are not the major allele of
/// its variant, so it is a whole number from 0 to the ploidy, and a
/// genotype with any allele missing takes the mean dosage of its variant
/// instead, which is 0 for a variant with no called genotype at all. The
/// major allele and that mean are of the tested individuals, as everything
/// else here is. The dosages are not divided by anything: `beta` of the
/// result is the effect of one more copy of a non major allele in the
/// units of the trait, which the deviation of the variant would turn into
/// deviations.
///
/// The buffers are made as long as a block needs and are kept from one
/// block to the next, so a pass over a million variants asks the machine
/// for them once and allocates nothing for a variant.
#[derive(Debug, Clone)]
pub(crate) struct GwasDosages {
    /// How many variants the block held, which is how many rows of the
    /// result it gives.
    pub(super) num_vars: usize,
    /// How many individuals the study tests, which is the length of one
    /// row of the dosages.
    num_individuals: usize,
    /// How many variants of the block have variance among them.
    pub(super) num_with_variance: usize,
    /// The dosages, with the rows of the variants that have variance at
    /// the start, in the order of the block.
    dosages: Vec<f64>,
    /// The frequency of each variant of the block, in its order.
    pub(super) allele_freq: Vec<f64>,
    /// Whether each variant of the block has variance, in its order.
    pub(super) has_variance: Vec<bool>,
}

impl GwasDosages {
    /// The buffers of a study that has read no block yet.
    #[must_use]
    pub(crate) fn of_a_study() -> GwasDosages {
        GwasDosages {
            num_vars: 0,
            num_individuals: 0,
            num_with_variance: 0,
            dosages: Vec::new(),
            allele_freq: Vec::new(),
            has_variance: Vec::new(),
        }
    }

    /// The dosages of the variants of `block` over the individuals the
    /// study tests, which replace whatever block was read before.
    ///
    /// The individuals that are not tested leave the block before anything
    /// is counted, so the major allele of a variant, the mean its missing
    /// genotypes take, its frequency and whether it has any variance are
    /// all of the tested individuals and of nobody else. That is what
    /// "What it gives" of `docs/specs/gwas.md` asks for, and it matters as
    /// soon as a phenotype leaves one individual out: the block is the
    /// whole panel and the study is of those that have a trait. The block
    /// is cut down in place, so the caller is left with the genotypes of
    /// the tested individuals and a block that is read twice is a block
    /// whose individuals are already gone.
    ///
    /// `design` carries those individuals, checked, and whether a variant
    /// with more than two alleles is read. `of_the_pass` is the shape the
    /// reader says its source has and where the block sits among the
    /// variants it has given.
    ///
    /// Every block is checked against the reader before a row of it is
    /// read, which is what `alleles_per_var_of` of [`crate::stats`] does
    /// for the two passes of that module, with the same errors: a pass
    /// reads the rows of every block as rows of one run over the variants,
    /// so a block of other individuals or of another ploidy is read one
    /// individual at the place of another, and a block of no variant is a
    /// reader that has stopped without saying so.
    ///
    /// # Errors
    ///
    /// Each of these but the last is a defect of the reader that gave the
    /// block or of the caller, and not of the dataset.
    /// [`Error::FieldsNotInTheBlock`] when the block holds no genotypes,
    /// [`Error::ReaderGaveABlockOfNoVariants`] when it holds no variant,
    /// [`Error::BlocksDoNotFitTogether`] when its individuals or its
    /// ploidy are not the reader's, [`Error::BlockWithNoGenotypeOfAVariant`]
    /// when it holds the genotypes of no individual, and
    /// [`Error::BlockArrayOfAnotherSize`] when its arrays are not of the
    /// size it states, which [`Block::check`] finds and which the rows
    /// that came out are counted against again.
    /// [`Error::BlockTooLarge`] when its individuals times its ploidy are
    /// more than a `usize` counts. What
    /// [`Block::retain_individuals`] refuses of the tested individuals,
    /// which is an individual the block has not and one that is there
    /// twice. And what reading the dosages of one row refuses, with
    /// [`Error::GwasVariantsTooLarge`] when the position of a variant is
    /// beyond what a `usize` counts.
    pub(crate) fn read_the_block(
        &mut self,
        block: &mut Block,
        design: &Design<'_>,
        of_the_pass: BlockOfThePass,
    ) -> Result<()> {
        let missing = Needs::GTS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        // The rows of a block are cut out of the sizes it states, so those
        // sizes are checked before anything is read.
        block.check()?;
        if block.num_vars == 0 {
            return Err(Error::ReaderGaveABlockOfNoVariants);
        }
        // The individuals and the ploidy are compared with the reader's
        // before the block is cut down to the tested individuals, since
        // that is what changes the first of the two. A block of another
        // ploidy is the one that no other check catches: its rows would be
        // cut at one width and its genotypes read at another, and when the
        // two disagree enough the rows come out as none at all.
        if block.num_individuals != design.num_individuals_of_the_source()
            || block.ploidy != of_the_pass.ploidy
        {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals: design.num_individuals_of_the_source(),
                ploidy: of_the_pass.ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        if block.alleles_per_var()? == 0 {
            return Err(Error::BlockWithNoGenotypeOfAVariant {
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
            });
        }
        if !design
            .individuals()
            .iter()
            .copied()
            .eq(0..block.num_individuals)
        {
            block.retain_individuals(design.individuals())?;
        }
        let (Some(num_individuals), Some(alleles_per_var)) = (
            NonZeroUsize::new(block.num_individuals),
            NonZeroUsize::new(block.alleles_per_var()?),
        ) else {
            // The block held the genotypes of one individual at least
            // above, and the tested individuals are one at least, since
            // `retain_individuals` refuses none: neither of these is 0.
            return Err(Error::BlockWithNoGenotypeOfAVariant {
                num_individuals: block.num_individuals,
                ploidy: block.ploidy,
            });
        };
        // The block is of the reader's ploidy and its genotypes are its
        // variants times the alleles of one of them, so this division is
        // exact and the buffer holds one value for each individual of each
        // row.
        let num_values = block
            .num_vars
            .checked_mul(num_individuals.get())
            .ok_or(Error::GwasVariantsTooLarge)?;
        self.dosages.resize(num_values, 0.0);
        let rows = the_dosages_of_the_rows(
            &block.gts,
            alleles_per_var,
            num_individuals,
            of_the_pass.ploidy,
            design.multiallelic,
            of_the_pass.first_var,
            &mut self.dosages,
        )?;
        // The rows are cut out of the genotypes at the width of one
        // variant, so a block whose genotypes are not its variants times
        // that width gives fewer rows than it says it holds, and the
        // variants that are left over would go out of the result with no
        // error. Everything above says that cannot happen here; this is
        // what says so of the rows that actually came out.
        if rows.len() != block.num_vars {
            return Err(Error::BlockArrayOfAnotherSize {
                array: "gts",
                found: block.gts.len(),
                expected: block.num_vars.saturating_mul(alleles_per_var.get()),
            });
        }
        self.num_vars = rows.len();
        self.num_individuals = num_individuals.get();
        self.allele_freq.clear();
        self.allele_freq
            .extend(rows.iter().map(|row| row.allele_freq));
        self.has_variance.clear();
        self.has_variance
            .extend(rows.iter().map(|row| row.has_variance));
        self.num_with_variance = self
            .has_variance
            .iter()
            .filter(|has_variance| **has_variance)
            .count();
        // The rows of the variants that have variance are moved to the
        // start of the buffer, so that the block is one matrix of the
        // variants a model can test and the ones that have no answer are
        // not in it. A block with no variant to leave out moves nothing.
        for (to, (var, _)) in self
            .has_variance
            .iter()
            .enumerate()
            .filter(|(_, has_variance)| **has_variance)
            .enumerate()
        {
            if to != var {
                let from = the_row_of(var, self.num_individuals);
                let start = the_row_of(to, self.num_individuals).start;
                self.dosages.copy_within(from, start);
            }
        }
        Ok(())
    }

    /// How many variants the block held, which is how many rows of the
    /// result it gives: the ones that have no answer are among them. It is
    /// 0 before a block has been read.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals the study tests, which is the length of one
    /// row of [`GwasDosages::dosages`]. It is 0 before a block has been
    /// read, since it is the block that says which of its individuals were
    /// kept.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many variants of the block have variance among the tested
    /// individuals, which is how many rows [`GwasDosages::dosages`] holds.
    #[must_use]
    pub(crate) fn num_with_variance(&self) -> usize {
        self.num_with_variance
    }

    /// The dosages of the variants that have variance, one row of
    /// [`GwasDosages::num_individuals`] values for each of them, in the
    /// order of the block. It is what a model tests as one matrix, and it
    /// is empty before a block has been read.
    ///
    /// A dosage is a whole number from 0 to the ploidy, how many alleles
    /// of the genotype are not the major allele of its variant among the
    /// tested individuals, and a genotype with any allele missing holds
    /// the mean dosage of its variant instead, which is what centering the
    /// variant would make 0.
    #[must_use]
    pub(crate) fn dosages(&self) -> &[f64] {
        // The buffer holds one row for every variant of the block, and the
        // rows of the variants that have variance were moved to its start,
        // so it holds this many values at least. A buffer that did not
        // would be a defect of this module, and what it gives then is no
        // value at all and not a longer slice, which a model would read as
        // more variants than the block holds.
        let values = self.num_with_variance.saturating_mul(self.num_individuals);
        self.dosages.get(..values).unwrap_or_default()
    }

    /// The frequency of the alleles that are not the major one, over the
    /// tested individuals: one for each variant of the block, in its
    /// order, the variants that have no answer among them. It is the mean
    /// dosage of the variant over the ploidy, and 0 for a variant with no
    /// called genotype among those individuals.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn allele_freq(&self) -> &[f64] {
        &self.allele_freq
    }

    /// Whether each variant of the block has variance among the tested
    /// individuals, in the order of the block. A variant that has none has
    /// no answer, as "The variants that have no answer" of
    /// `docs/specs/gwas.md` says.
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "what the block held, which the models that report it will read; \
                      the study itself calls `of_a_study`, `read_the_block`, \
                      `num_with_variance` and `dosages`, and until then the tests of \
                      this module are what read this"
        )
    )]
    pub(crate) fn has_variance(&self) -> &[bool] {
        &self.has_variance
    }
}

/// The dosages of a block over the individuals a study tests, which are
/// theirs and not the whole panel's, as "What it gives" of
/// `docs/specs/gwas.md` asks.
#[cfg(test)]
mod tests {
    use super::{BlockOfThePass, GwasDosages};
    use crate::error::Error;
    use crate::gwas::fixtures::{
        FREQUENCIES_OF_THE_PANEL, FREQUENCIES_OF_THE_TESTED, MISSING, OF_EIGHT, TESTED_OF_EIGHT,
        THE_PANEL_OF_EIGHT, a_block, a_study, assert_the_values, the_design_of, the_first_block_of,
        the_phenotype_and_the_design_of,
    };

    /// The dosages, the major allele, the mean a genotype with an allele
    /// missing takes and the frequency are all of the individuals that
    /// have a phenotype, and the four that have none change every one of
    /// them.
    ///
    /// The fixture is the panel of eight diploid individuals of
    /// `fixtures::OF_EIGHT`, of which the four at 0, 2, 4 and 6 are
    /// tested: the doc comments of `FREQUENCIES_OF_THE_TESTED` and of
    /// `FREQUENCIES_OF_THE_PANEL` work every frequency here out by hand
    /// from the genotypes, and this test asserts both, the study of the
    /// four and the study of all eight, so that neither is read off the
    /// other. Three of the four frequencies move when the individuals
    /// without a phenotype are read, by 0.25, by 0.048 and by 0.5, which
    /// are 2.5e14, 4.8e13 and 5e14 times the 1e-15 the values are held to.
    ///
    /// The fourth, `v0`, keeps its frequency and loses its answer instead:
    /// every tested individual is heterozygous, so it has no variance
    /// among them and cannot be tested, while over the panel it has
    /// variance and would be. A study that read the panel would give it
    /// three numbers where it has none.
    #[test]
    fn the_dosages_and_the_frequencies_are_of_the_tested_individuals_and_not_of_the_panel() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        // The four individuals that have a phenotype are four of the eight
        // the source has, in its order, so this is a study that runs, and
        // the design that comes out is what carries them to the block.
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);
        let mut dosages = GwasDosages::of_a_study();

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block over the four tested individuals");

        assert_eq!(dosages.num_vars(), 4);
        assert_eq!(dosages.num_individuals(), 4, "the tested individuals");
        assert_the_values(
            dosages.allele_freq(),
            &FREQUENCIES_OF_THE_TESTED,
            "the frequencies over the tested individuals",
        );
        assert_eq!(
            dosages.has_variance(),
            [false, true, true, false],
            "v0 is heterozygous in all four and v3 has no genotype among them"
        );
        assert_eq!(dosages.num_with_variance(), 2);
        // The rows of the two variants that have variance, in the order of
        // the block and at the start of the buffer: `v1` is 0/0 0/1 1/1
        // 0/1, whose dosages are 0 1 2 1, and `v2` is 0/0 ./. 1/1 0/0,
        // whose dosages are 0, the mean 2 / 3 of the other three, 2 and 0.
        assert_the_values(
            dosages.dosages(),
            &[0.0, 1.0, 2.0, 1.0, 0.0, 0.6666666666666666, 2.0, 0.0],
            "the dosages of the variants that have variance",
        );

        let (phenotype, design) = the_phenotype_and_the_design_of(&THE_PANEL_OF_EIGHT);
        let study = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &OF_EIGHT);

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the block over all eight individuals");

        assert_eq!(dosages.num_individuals(), 8, "the whole panel");
        assert_the_values(
            dosages.allele_freq(),
            &FREQUENCIES_OF_THE_PANEL,
            "the frequencies over the whole panel",
        );
        assert_eq!(
            dosages.has_variance(),
            [true, true, true, true],
            "every variant of the panel has variance, `v0` and `v3` among them"
        );
        // `v2` is where the two studies count the dosages from different
        // alleles: over the panel the major allele is 1, so the first
        // individual, 0/0, has the dosage 2, where over the tested four it
        // has 0.
        assert_the_values(
            dosages.dosages().get(16..24).unwrap_or_default(),
            &[2.0, 0.0, 0.5714285714285714, 0.0, 0.0, 0.0, 2.0, 0.0],
            "the dosages of `v2` over the whole panel",
        );
    }

    /// The three variants of "The worked example" of
    /// `docs/specs/gwas.md`, whose dosages and frequencies the spec gives,
    /// read over its six individuals while two more are in the block and
    /// have no phenotype.
    ///
    /// The spec's dosages are 0 1 2 0 1 2 for `v0`, 0 1 2 0.8 1 0 for
    /// `v1`, where the 0.8 is the genotype of `i3` that was not called
    /// taking the mean `(0 + 1 + 2 + 1 + 0) / 5` of its variant, and
    /// 1 1 1 1 1 1 for `v2`, which has no variance and is not kept; the
    /// frequencies are 0.5, 0.4 and 0.5.
    ///
    /// The two individuals that are left out are `1/1` at all three
    /// variants, which is not what the six hold: it makes 1 the major
    /// allele of every one of the three, so every dosage is counted from
    /// the other allele, and it gives `v2` two dosages and an answer it
    /// does not have among the six. Over all eight the frequencies are
    /// 0.375, 3 / 7 and 0.375, so a study that read them misses all three
    /// of the spec's numbers by 0.125, 0.029 and 0.125.
    #[test]
    fn the_dosages_of_the_worked_example_are_of_its_six_individuals_and_not_of_the_two_beside_them()
    {
        let of_the_worked_example: [&[i8]; 3] = [
            // v0: 0/0 0/1 1/1 0/0 0/1 1/1 | 1/1 1/1
            &[0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
            // v1: 0/0 0/1 1/1 ./. 0/1 0/0 | 1/1 1/1
            &[0, 0, 0, 1, 1, 1, MISSING, MISSING, 0, 1, 0, 0, 1, 1, 1, 1],
            // v2: 0/1 0/1 0/1 0/1 0/1 0/1 | 1/1 1/1
            &[0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1],
        ];
        let six = [0, 1, 2, 3, 4, 5];
        let (phenotype, design) = the_phenotype_and_the_design_of(&six);
        let study = a_study(&phenotype, &design, &six);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &of_the_worked_example);
        let mut dosages = GwasDosages::of_a_study();

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the worked example");

        assert_the_values(
            dosages.allele_freq(),
            &[0.5, 0.4, 0.5],
            "the frequencies of the worked example",
        );
        assert_eq!(dosages.has_variance(), [true, true, false]);
        assert_the_values(
            dosages.dosages(),
            &[0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.8, 1.0, 0.0],
            "the dosages of the worked example",
        );

        let eight = [0, 1, 2, 3, 4, 5, 6, 7];
        let (phenotype, design) = the_phenotype_and_the_design_of(&eight);
        let study = a_study(&phenotype, &design, &eight);
        let design = the_design_of(&study, 8);
        let mut block = a_block(8, 2, &of_the_worked_example);

        dosages
            .read_the_block(&mut block, &design, the_first_block_of(2))
            .expect("the dosages of the worked example and the two beside it");

        assert_the_values(
            dosages.allele_freq(),
            &[0.375, 0.42857142857142855, 0.375],
            "the frequencies of the eight individuals of the block",
        );
        assert_eq!(
            dosages.has_variance(),
            [true, true, true],
            "`v2` has variance once the two individuals that are 1/1 are read"
        );
    }

    /// How many alleles a variant has is counted among the tested
    /// individuals too, so a third allele that only an individual without
    /// a phenotype carries is not a variant of more than two alleles for
    /// the study, and one a tested individual carries is refused with its
    /// position among the variants the reader has given.
    ///
    /// The fixture is six triploid individuals, of which the four at 1, 2,
    /// 3 and 5 are tested, and two variants, the second of which has the
    /// allele 2 in the individual 0, who has no phenotype. The position
    /// asserted is 101 and not 1, because the block is given as the
    /// hundred and first variant of the reader.
    ///
    /// With `transform_to_biallelic` the whole panel is read and every
    /// allele that is not the major one counts the same: the alleles of
    /// the second variant over the six are fourteen 0s, three 1s and one
    /// 2, so the major one is 0 and the dosages are 1 0 0 1 0 2, whose
    /// mean is 4 / 6 and whose frequency over the ploidy of 3 is
    /// 2 / 9 = 0.2222222222222222.
    #[test]
    fn a_third_allele_is_counted_among_the_tested_individuals_and_refused_there() {
        let of_six_triploids: [&[i8]; 2] = [
            // w0: 0/0/0 0/0/1 1/1/1 0/1/1 0/0/0 1/1/1
            &[0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1],
            // w1: 0/0/2 0/0/0 0/0/0 0/0/1 0/0/0 0/1/1
            &[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1],
        ];
        let four = [1, 2, 3, 5];
        let (phenotype, design) = the_phenotype_and_the_design_of(&four);
        let study = a_study(&phenotype, &design, &four);
        let of_four = the_design_of(&study, 6);
        let mut block = a_block(6, 3, &of_six_triploids);
        let mut dosages = GwasDosages::of_a_study();
        let of_the_pass = BlockOfThePass {
            ploidy: 3,
            first_var: 100,
        };

        dosages
            .read_the_block(&mut block, &of_four, of_the_pass)
            .expect("the allele 2 is the individual 0's, who has no phenotype");

        // `w0` over the four tested individuals is 0/0/1 1/1/1 0/1/1
        // 1/1/1, three 0s and nine 1s, so the major allele is 1 and the
        // dosages are 2 0 1 0, whose mean is 3 / 4 and whose frequency
        // over the ploidy of 3 is 0.25. `w1` is 0/0/0 0/0/0 0/0/1 0/1/1,
        // nine 0s and three 1s, the major allele 0 and the dosages
        // 0 0 1 2, the same mean and the same frequency.
        assert_the_values(
            dosages.allele_freq(),
            &[0.25, 0.25],
            "the frequencies over the four tested individuals",
        );
        assert_the_values(
            dosages.dosages(),
            &[2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2.0],
            "the dosages over the four tested individuals",
        );

        let six = [0, 1, 2, 3, 4, 5];
        let (phenotype, values) = the_phenotype_and_the_design_of(&six);
        let study = a_study(&phenotype, &values, &six);
        let of_six = the_design_of(&study, 6);
        let mut block = a_block(6, 3, &of_six_triploids);

        match dosages.read_the_block(&mut block, &of_six, of_the_pass) {
            Err(Error::VariantWithMoreThanTwoAlleles {
                position,
                num_alleles,
            }) => {
                assert_eq!((position, num_alleles), (101, 3));
            }
            other => panic!(
                "the individual 0 carries the allele 2 of the second variant of the \
                 block, which is the variant 101 of the reader, and that gave {other:?}"
            ),
        }

        let mut collapsed = a_study(&phenotype, &values, &six);
        collapsed.transform_to_biallelic = true;
        let collapsed = the_design_of(&collapsed, 6);
        let mut block = a_block(6, 3, &of_six_triploids);

        dosages
            .read_the_block(&mut block, &collapsed, of_the_pass)
            .expect("every allele that is not the major one counts the same");

        assert_the_values(
            dosages.allele_freq(),
            &[0.5, 0.2222222222222222],
            "the frequencies over the six individuals with the third allele collapsed",
        );
    }

    /// The rows of a block are read on the threads of rayon, and neither
    /// the dosages nor the frequencies depend on how many there are, nor
    /// on whether the rows were read one after another, which is what
    /// WebAssembly does.
    ///
    /// The block is 300 variants of 40 diploid individuals, of which the
    /// 14 whose position is a multiple of 3 are tested, so a pool of four
    /// shares the rows out in several chunks. The alleles of the genotype
    /// of the individual `k` at the variant `v` are read out of a pattern
    /// of 5 at `(7k + 3v) % 5` and one of 7 at `(11k + 5v) % 7`, both of
    /// them the missing allele, 0 and 1, and every fiftieth variant is
    /// `0/0` in everybody instead.
    ///
    /// What that holds, counted over the 14 tested individuals on 23
    /// September 2026: they take 6 to 8 different genotypes at each of the
    /// variants that are not the fiftieth, of the 9 the two patterns can
    /// give, and never fewer, since a tested individual is a multiple of 3
    /// and the two patterns then run over all of 5 and all of 7; 30.8 per
    /// cent of their genotypes have an allele missing and no variant is
    /// missing in all of them; no row is the row of its neighbour, and the
    /// 300 rows are 36 different ones, so a row read at the place of
    /// another is seen; and the 6 variants that are `0/0` in everybody
    /// have no variance among them, so the two runs have to agree about
    /// which rows are left out of the matrix and where the others moved
    /// to. The fixture before this one gave the tested individuals 3
    /// genotypes and 5 of them the same one at every variant.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those
    /// alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_rows_are_the_same_on_one_thread_on_several_and_read_one_after_another() {
        use std::num::NonZeroUsize;

        use super::{
            MultiallelicVariants, the_dosages_of_the_rows, the_dosages_of_the_rows_one_by_one,
        };

        /// Whether two lists of values hold the same numbers, to the bit,
        /// which is what two ways of reading the same genotypes give: each
        /// row is read on its own and nothing is summed across the rows,
        /// so neither the threads nor their number can move a digit.
        fn the_same_values(one: &[f64], other: &[f64]) -> bool {
            one.len() == other.len()
                && one
                    .iter()
                    .zip(other)
                    .all(|(one, other)| one.to_bits() == other.to_bits())
        }

        // The missing allele, 0 and 1 in a pattern of 5 and one of 7,
        // which no two of the tested individuals read at the same place.
        const OF_THE_FIRST: [i8; 5] = [0, 1, MISSING, 1, 0];
        const OF_THE_SECOND: [i8; 7] = [0, 1, 0, 1, MISSING, 1, 0];
        let genotype = |individual: usize, variant: usize| match variant % 50 {
            0 => [0, 0],
            _ => [
                OF_THE_FIRST[(individual * 7 + variant * 3) % 5],
                OF_THE_SECOND[(individual * 11 + variant * 5) % 7],
            ],
        };
        let of_forty = || {
            let rows: Vec<Vec<i8>> = (0..300)
                .map(|variant| {
                    (0..40)
                        .flat_map(|individual| genotype(individual, variant))
                        .collect()
                })
                .collect();
            let rows: Vec<&[i8]> = rows.iter().map(Vec::as_slice).collect();
            a_block(40, 2, &rows)
        };
        let tested: Vec<usize> = (0..40).filter(|individual| individual % 3 == 0).collect();
        let (phenotype, design) = the_phenotype_and_the_design_of(&tested);
        let study = a_study(&phenotype, &design, &tested);
        let design = the_design_of(&study, 40);
        let read_with = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            let mut block = of_forty();
            let mut dosages = GwasDosages::of_a_study();
            pool.install(|| dosages.read_the_block(&mut block, &design, the_first_block_of(2)))
                .expect("the dosages of the block");
            dosages
        };

        let on_one = read_with(1);
        let on_four = read_with(4);

        assert_eq!(on_one.num_vars(), 300);
        assert_eq!(on_one.num_individuals(), 14);
        assert_eq!(
            on_one.num_with_variance(),
            294,
            "the 6 variants that are 0/0 in everybody have no variance"
        );
        assert!(
            the_same_values(on_one.dosages(), on_four.dosages()),
            "the dosages on one thread and on four"
        );
        assert!(
            the_same_values(on_one.allele_freq(), on_four.allele_freq()),
            "the frequencies on one thread and on four"
        );
        assert_eq!(on_one.has_variance(), on_four.has_variance());

        // The rows read one after another, which is the pass WebAssembly
        // takes, against the rows read on the threads, over the genotypes
        // of the tested individuals alone: the block is compacted to them
        // here as the two passes above had it compacted for them.
        let mut block = of_forty();
        block
            .retain_individuals(&tested)
            .expect("the tested individuals");
        let alleles_per_var =
            NonZeroUsize::new(block.alleles_per_var().expect("the alleles of a variant"))
                .expect("the block holds genotypes");
        let num_individuals = NonZeroUsize::new(block.num_individuals).expect("the individuals");
        let mut on_threads = vec![0.0; block.gts.len() / 2];
        let rows = the_dosages_of_the_rows(
            &block.gts,
            alleles_per_var,
            num_individuals,
            2,
            MultiallelicVariants::Refused,
            0,
            &mut on_threads,
        )
        .expect("the rows on the threads");
        let mut one_by_one = vec![0.0; block.gts.len() / 2];
        let rows_one_by_one = the_dosages_of_the_rows_one_by_one(
            &block.gts,
            alleles_per_var,
            num_individuals,
            2,
            MultiallelicVariants::Refused,
            0,
            &mut one_by_one,
        )
        .expect("the rows one after another");

        assert!(
            the_same_values(&on_threads, &one_by_one),
            "the dosages of every row, on the threads and one after another"
        );
        assert_eq!(rows.len(), 300);
        for (var, (row, one_by_one)) in rows.iter().zip(&rows_one_by_one).enumerate() {
            assert_eq!(
                (row.allele_freq.to_bits(), row.has_variance),
                (one_by_one.allele_freq.to_bits(), one_by_one.has_variance),
                "the variant {var}"
            );
        }
    }

    /// The dosages of a study and the standardized row of
    /// [`crate::variant`] read a variant by the same three rules, and this
    /// is what fails the day one of them is changed and the other is not.
    ///
    /// The two are written separately because a study needs what the other
    /// one does not give: the mean of the variant, which it reports as
    /// `allele_freq`, and dosages that are not divided by the deviation of
    /// the variant, since `beta` is the effect of one more copy of an
    /// allele in the units of the trait. Everything before that division
    /// is the same work, and this asserts that it gives the same answers:
    /// which allele is the major one, and so what every dosage is counted
    /// from; the mean a genotype with an allele missing takes; whether the
    /// variant has variance at all; and the refusal of a variant with more
    /// than two alleles, at the same position and with the same count.
    ///
    /// The fixture is the panel of eight over its four tested individuals,
    /// which holds a variant with no variance among them, two that have
    /// it, one of which has a genotype that was not called, and one with
    /// nothing called at all. Each row of it is read both ways, and the
    /// standardized row is asserted to be this one centered at the mean
    /// this module gives and divided by the deviation worked out from it,
    /// which is `docs/specs/pca.md`'s divisor: the root of the mean square
    /// deviation over all the individuals, the ones whose genotype was not
    /// called counting as no deviation, since they hold the mean.
    ///
    /// The bound is 1e-12 of the largest value of the row, which is the
    /// scale of what is being compared, and not of each value, which is 0
    /// for an individual at the mean. The two paths add the same squares
    /// in a different order, this one over the individuals and the other
    /// over the dosages with a count on each, so the last bits may differ;
    /// measured over this fixture on 23 September 2026 the worst distance
    /// was 0 on Accelerate and on faer, with the bound set to 0.
    #[test]
    fn the_dosages_and_the_standardized_row_of_a_variant_agree() {
        use super::{DosageScratch, MultiallelicVariants, the_dosages_of_a_row};
        use crate::variant::{DosageOptions, DosageScale, RowScratch, the_standardized_row};

        let options = DosageOptions {
            transform_to_biallelic: false,
            scale: DosageScale::OfTheDosages,
        };
        let mut block = a_block(8, 2, &OF_EIGHT);
        block
            .retain_individuals(&TESTED_OF_EIGHT)
            .expect("the four tested individuals");
        let mut mine = DosageScratch::of(4);
        let mut theirs = RowScratch::of(4);
        let mut agreed: Vec<bool> = Vec::new();

        for (var, gts) in block.gts.as_chunks::<8>().0.iter().enumerate() {
            let mut row = [0.0; 4];
            let read = the_dosages_of_a_row(
                gts,
                2,
                var,
                MultiallelicVariants::Refused,
                &mut mine,
                &mut row,
            )
            .expect("the dosages of the row");
            let mut standardized = [0.0; 4];
            let used = the_standardized_row(gts, 2, var, &options, &mut theirs, &mut standardized)
                .expect("the standardized row");

            assert_eq!(
                read.has_variance, used,
                "the variant {var} has variance one way and not the other"
            );
            agreed.push(used);
            if !used {
                continue;
            }
            // The mean of the called dosages, which is the frequency times
            // the ploidy: the two multiply and divide by 2, which is exact.
            let mean = read.allele_freq * 2.0;
            let squares: f64 = row
                .iter()
                .map(|dosage| (dosage - mean) * (dosage - mean))
                .sum();
            let divisor = (squares / 4.0).sqrt();
            let largest = standardized
                .iter()
                .fold(0.0_f64, |largest, value| largest.max(value.abs()));
            for (position, (theirs, mine)) in standardized.iter().zip(&row).enumerate() {
                let from_the_dosage = (mine - mean) / divisor;
                assert!(
                    (from_the_dosage - theirs).abs() <= 1e-12 * largest,
                    "the variant {var} of the individual {position}: the dosage {mine} \
                     centered and divided is {from_the_dosage} and the standardized row \
                     holds {theirs}"
                );
            }
        }

        assert_eq!(
            agreed,
            [false, true, true, false],
            "the variants of the panel that have variance among the four tested"
        );

        // A variant of three alleles: one genotype of the individual 0 is
        // 0/0/2 and the others hold 0 and 1, and both refuse it with its
        // position among the variants the reader gave and the count.
        let block = a_block(
            6,
            3,
            &[&[0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1]],
        );
        let mut mine = DosageScratch::of(6);
        let mut theirs = RowScratch::of(6);
        let mut row = [0.0; 6];
        let mut standardized = [0.0; 6];

        let refused = the_dosages_of_a_row(
            &block.gts,
            3,
            7,
            MultiallelicVariants::Refused,
            &mut mine,
            &mut row,
        );
        let refused_there =
            the_standardized_row(&block.gts, 3, 7, &options, &mut theirs, &mut standardized);

        match (refused, refused_there) {
            (
                Err(Error::VariantWithMoreThanTwoAlleles {
                    position,
                    num_alleles,
                }),
                Err(Error::VariantWithMoreThanTwoAlleles {
                    position: there,
                    num_alleles: alleles_there,
                }),
            ) => {
                assert_eq!((position, num_alleles), (7, 3));
                assert_eq!((there, alleles_there), (7, 3));
            }
            (mine, theirs) => panic!(
                "the variant has three alleles, and the dosages gave {mine:?} and the \
                 standardized row gave {theirs:?}"
            ),
        }
    }

    /// A block that holds no genotype is refused, and so is a block of no
    /// variant, which is a reader that has stopped without saying so.
    #[test]
    fn a_block_that_holds_no_genotype_and_a_block_of_no_variant_are_refused() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);
        block.gts = Vec::new();

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::FieldsNotInTheBlock { fields }) => {
                assert_eq!(fields, crate::variant::Needs::GTS);
            }
            other => panic!("the block holds no genotype, and that gave {other:?}"),
        }

        let mut block = a_block(8, 2, &[]);

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::ReaderGaveABlockOfNoVariants) => {}
            other => panic!("the block holds no variant, and that gave {other:?}"),
        }
    }

    /// A block whose arrays are not of the size it states is refused when
    /// every individual is tested, which is when nothing else looks at it:
    /// the rows are cut out of the genotypes by that size, so a block that
    /// is one allele short is read as three variants where it says four,
    /// and the fourth would leave the result with no error.
    ///
    /// A study of some of the individuals reaches the same refusal through
    /// `Block::retain_individuals`, which checks the block itself before
    /// it moves an allele. A study of all of them cuts nothing down, so
    /// this is the case that says the check is made here as well.
    #[test]
    fn a_block_whose_arrays_are_not_of_its_size_is_refused_with_every_individual_tested() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&THE_PANEL_OF_EIGHT);
        let study = a_study(&phenotype, &design, &THE_PANEL_OF_EIGHT);
        let design = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);
        block.gts.pop();

        match dosages.read_the_block(&mut block, &design, the_first_block_of(2)) {
            Err(Error::BlockArrayOfAnotherSize {
                array,
                found,
                expected,
            }) => {
                assert_eq!((array, found, expected), ("gts", 63, 64));
            }
            other => panic!("the block is one allele short, and that gave {other:?}"),
        }
    }

    /// A block of another ploidy or of other individuals than the reader
    /// says its source has is refused, naming both shapes.
    ///
    /// The ploidy is the one no other check catches. The rows of a block
    /// are cut at its own individuals times its own ploidy, and the buffer
    /// they are read into is sized at the pass's ploidy, so when the two
    /// disagree the rows and the buffer do not line up. Both of these were
    /// run against the commit before this one on 23 September 2026: the
    /// haploid variant of five individuals below, read at the ploidy 2,
    /// gave `Ok` with every variant of the block gone and nothing to show
    /// it, because the buffer came out shorter than one row and the two
    /// zipped to nothing, so no row ran at all; the diploid block of eight
    /// read at the ploidy 6 got as far as cutting a row and stopped there,
    /// with an error about genotypes that are not whole, which names the
    /// genotypes of the block and not the ploidy the pass was reading at.
    /// Which of the two a block gets turns on how far the sizes are apart
    /// and on how many individuals are tested, and neither is an answer.
    /// The frequency of a variant is its mean dosage over the pass's
    /// ploidy, so even a block that came out whole would be read into
    /// frequencies of nothing.
    #[test]
    fn a_block_of_another_ploidy_or_of_other_individuals_than_the_readers_is_refused() {
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let of_eight = the_design_of(&study, 8);
        let mut dosages = GwasDosages::of_a_study();
        let mut block = a_block(8, 2, &OF_EIGHT);

        match dosages.read_the_block(&mut block, &of_eight, the_first_block_of(6)) {
            Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                ploidy,
                found_num_individuals,
                found_ploidy,
            }) => {
                assert_eq!((num_individuals, ploidy), (8, 6));
                assert_eq!((found_num_individuals, found_ploidy), (8, 2));
            }
            other => panic!("the block is diploid and the pass reads 6, and that gave {other:?}"),
        }

        // One haploid variant of five individuals, four of which are
        // tested, read by a pass that says its source is diploid.
        let four = [0, 1, 3, 4];
        let (phenotype, design) = the_phenotype_and_the_design_of(&four);
        let study = a_study(&phenotype, &design, &four);
        let of_five = the_design_of(&study, 5);
        let mut block = a_block(5, 1, &[&[0, 1, 1, 1, 0]]);

        match dosages.read_the_block(&mut block, &of_five, the_first_block_of(2)) {
            Err(Error::BlocksDoNotFitTogether {
                ploidy,
                found_ploidy,
                ..
            }) => {
                assert_eq!((ploidy, found_ploidy), (2, 1));
            }
            other => panic!("the block is haploid and the pass reads 2, and that gave {other:?}"),
        }

        // The same four individuals of a source the reader says has nine,
        // which is a block of others: the positions of the tested
        // individuals are positions among the reader's.
        let (phenotype, design) = the_phenotype_and_the_design_of(&TESTED_OF_EIGHT);
        let study = a_study(&phenotype, &design, &TESTED_OF_EIGHT);
        let of_nine = the_design_of(&study, 9);
        let mut block = a_block(8, 2, &OF_EIGHT);

        match dosages.read_the_block(&mut block, &of_nine, the_first_block_of(2)) {
            Err(Error::BlocksDoNotFitTogether {
                num_individuals,
                found_num_individuals,
                ..
            }) => {
                assert_eq!((num_individuals, found_num_individuals), (9, 8));
            }
            other => panic!(
                "the reader says its source has nine individuals and the block has \
                 eight, and that gave {other:?}"
            ),
        }
    }
}
