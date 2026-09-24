//! Linkage disequilibrium: how much the genotype of one variant says about
//! the genotype of another.
//!
//! Two variants are in linkage disequilibrium when they sit close enough on
//! a chromosome that few recombinations have separated them, and the
//! measure of it is r², the square of the correlation between the dosages
//! of the two variants over the individuals called at both.
//! `docs/specs/ld.md` is the spec of the module.
//!
//! The module gives r² in the two shapes that spec asks for.
//! [`calc_r2_matrix`] reads a pass to its end and gives the r² of every
//! pair of its variants as a square matrix, with the chromosome and the
//! position of each variant. It is what a plot of one region is drawn
//! from, and how many variants it takes is bounded by the memory of the
//! matrix, which holds one value for each pair and so grows with the
//! square of them. [`calc_ld_and_dist`] gives how r² falls off as the two
//! variants of a pair move apart along a chromosome, for each population
//! of the dataset on its own, in bins of distance. It holds only the
//! variants that are still within reach of the newest one read, so it runs
//! over a dataset whose matrix of every pair no machine would hold.
//!
//! The dosage of a genotype is how many of its alleles are not the major
//! allele of its variant, 0, 1 or 2 in a diploid, and a genotype with an
//! allele missing has none. [`LdDosages`] reads the genotypes of a block as
//! the three matrices of variants x individuals that the six products of r²
//! are taken over: the dosages with a missing genotype written as 0, a 1
//! where the genotype was called and a 0 where it was not, and the square
//! of each dosage. Everything this module computes is taken from those
//! three, so the genotypes are read once for a whole set of variants and
//! never again.
//!
//! [`r2_between`] gives the r² of every variant of one set of dosages
//! against every variant of another, which is six products of those
//! matrices and the formula of the spec over the six sums they give.

use std::fmt;
use std::num::NonZeroUsize;

use popnei_linalg::{TheFirstOperand, TheSecondOperand, product};

use crate::block::Block;
use crate::block::BlockReader;
use crate::error::{Error, Result};
use crate::variant::{
    AlleleCounts, ChromTable, MISSING_ALLELE, Needs, count_alleles, the_major_allele,
    the_major_allele_frequency,
};

mod dist;

pub use dist::{
    DEFAULT_MAX_ALLOWED_MAF, DEFAULT_MAX_DIST, DEFAULT_MIN_DIST, DEFAULT_NUM_DIST_BINS, LdAndDist,
    LdAndDistOptions, LdBins, calc_ld_and_dist,
};

/// The most values one of the matrices of [`LdDosages`] holds, the variants
/// it was built over times the individuals.
///
/// It is what the routines of BLAS and LAPACK count the values of a matrix
/// in, and `crates/popnei-linalg` refuses a matrix above it on both of its
/// backends, so dosages that no product could be taken over are refused
/// where they are built and not at the first product. It is that crate's
/// [`THE_MOST_VALUES_OF_A_MATRIX`](popnei_linalg::THE_MOST_VALUES_OF_A_MATRIX)
/// under the name this module reads it by.
pub const MAX_VALUES_OF_THE_DOSAGES: usize = popnei_linalg::THE_MOST_VALUES_OF_A_MATRIX;

/// The most alleles a genotype of [`LdDosages`] holds.
///
/// A dosage is how many alleles of a genotype are not the major allele of
/// its variant, so it is at most the ploidy, and [`LdDosages::dosages`]
/// gives it as one byte. The VCF reader takes 255 alleles in a genotype at
/// most, and the largest ploidy of an organism is a dozen, so this refuses
/// no dataset that a reader of popnei gives.
pub const MAX_PLOIDY_OF_THE_DOSAGES: usize = 255;

/// The most alleles one variant of a set of dosages holds, its individuals
/// times the ploidy, for the r² to carry the bits the formula gives.
///
/// The six sums of a pair are whole numbers, and an `f64` holds each of
/// them exactly, but what the r² needs is that the four products the
/// formula takes of them be exact too: n·Σxx, n·Σxy, Σx·Σy and (Σx)². The
/// largest of the four is at most N²k², for N individuals of the ploidy k,
/// and a product of two whole numbers is exact while it is at most the
/// 2^53 up to which an `f64` counts one by one. So the bound is
/// Nk ≤ sqrt(2^53), which is this number: 47453132 diploid individuals,
/// and 372181 at the ploidy of 255 that the dosages take. Above it an r²
/// loses digits with nothing to show for it, 1.3e-12 of relative error at
/// a million individuals of the ploidy 255, which is past the 1e-12 the
/// spec compares within.
///
/// No dataset of this world reaches it: the objectives of popnei go to
/// 10000 individuals, and the largest ploidy of an organism is a dozen.
pub const MAX_ALLELES_OF_A_VARIANT: usize = 94_906_265;

/// How many variants [`calc_r2_matrix`] takes when the user names no
/// number, which is the default of `calc_rogers_huff_r2_matrix` in Python
/// and of `calcRogersHuffR2Matrix` in TypeScript.
///
/// The matrix holds one r² for each pair of the variants of the pass, so
/// it grows with the square of them: 200 MB of `f64` at this number and 80
/// GB at 100000. It is the one result of popnei that grows with the square
/// of its input, against goal 5 of `docs/objectives.md`, under which a
/// dataset never has to fit in memory, so a pass of more variants is
/// refused instead of asking the machine for the matrix of them. "Its
/// Python function" of `docs/specs/ld.md` gives the number; pyNei has
/// none and builds the matrix of whatever it is given.
pub const MAX_NUM_VARS_OF_THE_MATRIX: usize = 5000;

/// How many variants of the matrix one tile of the products holds.
///
/// The matrix is taken tile pair by tile pair and not over the whole set
/// at once, because each of the six sums of "How it runs" of
/// `docs/specs/ld.md` is a matrix over the pairs of the two tiles: six of
/// the whole set would be six times the result, 1.2 GB at 5000 variants,
/// where six of a pair of tiles of this size are 3 MB.
///
/// It was measured on 23 September 2026 by the performance review of
/// `docs/reports/perf-ld-2026-09-23.md`, on the owner's Apple M5 Pro with
/// the products on Accelerate, over the matrix of 5000 variants of 1000
/// individuals of the 400 MB VCF of `docs/rust_core.md`, the best of 5
/// runs of `crates/popnei/benches/r2_matrix.rs` with the reading of the
/// file taken out: 0.539 s at 128 variants, 0.449 s at 256, 0.424 s at
/// 512, 0.382 s at 1000, 0.380 s at 1250, 0.403 s at 2500 and 0.427 s at
/// 5000, which is one tile. Two things pull against each other. A larger
/// tile gives Accelerate a larger product, which it works out faster per
/// pair of variants, and it calls the linear algebra crate fewer times,
/// which scans both operands of every call for a value that is not
/// finite. Against that, a tile against itself computes its whole square
/// where the matrix needs half of it, so a larger tile computes more pairs
/// it throws away: 12819520 of them at 128 and 25000000 at 5000. The two
/// meet in a broad flat bottom from 1000 to 1250.
///
/// 1000 rather than 1250, which is 0.002 s faster at 5000 variants and
/// inside the spread of the runs, because the bottom moves with the
/// variants of the matrix and 1000 is the better of the two away from the
/// cap: over 4000 variants 1000 takes 0.255 s against 0.263 s at 1250, and
/// over 3000 it takes 0.147 s against 0.150 s. What those three sizes have
/// in common is that 1000 divides all of them, and a tile that does not
/// divide the variants leaves a short last tile whose products are shaped
/// badly: 1024 takes 0.414 s over 5000 variants against the 0.382 s of
/// 1000, 8 per 100 slower for a tile 2 per 100 larger, and it is the only
/// one of the sizes measured that leaves a remainder there.
///
/// The memory this costs is the six sums of one pair of tiles, which are
/// 6 x 8 bytes for each pair of variants of the pair: 3 MB at 256 and
/// 48 MB at 1000, beside the 200 MB of the matrix and the 120 MB of the
/// three matrices of the variants. In WebAssembly, where a `usize` is 32
/// bits and a tab holds a few GB, that is the cost to weigh if this
/// number is ever raised further.
///
/// The matrix does not change with it: the tiles cut the variants and
/// every sum of a pair runs over the individuals, so every sum is the same
/// whole number in whichever tile it is worked out, and
/// `neither_the_blocks_nor_the_tiles_change_the_matrix` compares the
/// matrix to the bit at tiles of 7, 64, 256, 500 and whatever this
/// constant says, against the matrix of one block.
const THE_VARS_OF_A_TILE: usize = 1000;

/// How the individuals of two sets of dosages whose r² was asked for
/// differ.
///
/// The sums of a pair run over the individuals both of its variants were
/// called in, which is the value at the same place of the two sets, so the
/// r² of two sets is taken only when both were built over the same
/// individuals of the block in the same order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TheIndividualsThatDiffer {
    /// The two sets were built over a different number of individuals.
    NotAsMany {
        /// How many individuals the first set was built over.
        of_a: usize,
        /// How many the second was built over.
        of_b: usize,
    },
    /// The two sets hold a different individual of the block at the same
    /// place.
    NotTheSame {
        /// The first place at which they differ, counted from 0.
        at: usize,
        /// The individual the first set holds there.
        of_a: usize,
        /// The individual the second set holds there.
        of_b: usize,
    },
}

impl fmt::Display for TheIndividualsThatDiffer {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::NotAsMany { of_a, of_b } => write!(
                formatter,
                "one was built over {of_a} individuals and the other over {of_b}"
            ),
            Self::NotTheSame { at, of_a, of_b } => write!(
                formatter,
                "the individual {at} of one is the individual {of_a} of the block and that of the other is the individual {of_b}"
            ),
        }
    }
}

/// The dosages of the variants of a block, for the individuals of one
/// population, held as the three matrices that the products of r² read.
///
/// Each matrix holds one value for every variant and individual, the
/// variants one after another, and the three are the A, the M and the S of
/// "How it runs" of `docs/specs/ld.md`: the dosages with a missing genotype
/// written as 0, a 1 where the genotype was called and a 0 where it was
/// not, and the square of each entry of the first. The six sums that r²
/// needs are the products of two of them, so the whole matrix of a set of
/// variants is six products and not one pass for each pair.
#[derive(Debug)]
pub struct LdDosages {
    /// How many variants the dosages are of.
    num_vars: usize,
    /// The individuals each variant has a value for, in the order they
    /// were given: the index each one has among the individuals of the
    /// block. [`LdDosages::of_block`] over an empty slice holds every
    /// individual of the block here, in the order the block has them, so
    /// that set and one built over all of them by name are the same set.
    individuals: Vec<usize>,
    /// The dosage of each genotype, with a genotype that has an allele
    /// missing written as 0. It is the A of the spec, and every entry of
    /// it is a whole number from 0 to the ploidy, which the dosages and
    /// the products of r² are worked out from.
    dosages: Vec<f64>,
    /// A 1 where the genotype was called and a 0 where an allele of it was
    /// not. It is the M of the spec, and the products over it count the
    /// individuals that a pair of variants has.
    called: Vec<f64>,
    /// The square of each entry of `dosages`. It is the S of the spec, and
    /// the products over it give the sums of the squares of a pair.
    squares: Vec<f64>,
    /// Whether the called genotypes of each variant hold two dosages at
    /// least.
    has_variance: Vec<bool>,
    /// The major allele frequency of each variant over these individuals,
    /// and `None` for a variant with no called allele.
    maf: Vec<Option<f64>>,
}

impl LdDosages {
    /// The dosages of the variants of `block`.
    ///
    /// `individuals` are indices into the individuals of the block, in the
    /// order they are given and each of them once, and an empty slice is
    /// every individual of it. The major allele of each variant is that of
    /// [`the_major_allele`](crate::variant::the_major_allele) over those
    /// individuals alone, so the dosages of a population are counted from
    /// the allele that population was called most often at.
    ///
    /// # Errors
    ///
    /// A block that does not pass [`Block::check`], one with variants and
    /// no genotypes, an index that is not an individual of the block, an
    /// individual asked for more than once, what the counts of one variant
    /// refuse, which are a variant of more alleles than a count of them
    /// holds and an allele below the missing one, a block whose
    /// individuals times its ploidy are more than
    /// [`MAX_ALLELES_OF_A_VARIANT`], a matrix this machine has not the
    /// memory for, a block whose genotypes hold more than
    /// [`MAX_PLOIDY_OF_THE_DOSAGES`] alleles each, and a block whose
    /// variants times its individuals is more than
    /// [`MAX_VALUES_OF_THE_DOSAGES`], which is what the linear algebra
    /// counts the values of a matrix in.
    pub fn of_block(block: &Block, individuals: &[usize]) -> Result<LdDosages> {
        block.check()?;
        let missing = Needs::GTS.difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        if block.ploidy > MAX_PLOIDY_OF_THE_DOSAGES {
            return Err(Error::LdPloidyTooLarge {
                ploidy: block.ploidy,
            });
        }
        let num_individuals = match individuals.is_empty() {
            true => block.num_individuals,
            false => individuals.len(),
        };
        // The alleles of one variant are counted before anything is
        // allocated: what this refuses is more individuals than any
        // machine holds the dosages of anyway.
        if num_individuals
            .checked_mul(block.ploidy)
            .is_none_or(|alleles| alleles > MAX_ALLELES_OF_A_VARIANT)
        {
            return Err(Error::LdTooManyAllelesInAVariant {
                num_individuals,
                ploidy: block.ploidy,
            });
        }
        // Which individuals have been asked for already, for the call
        // that names them: a call over every individual of the block names
        // none, and none can be named twice.
        let mut asked_for_already = match individuals.is_empty() {
            true => Vec::new(),
            false => a_vector_of(false, block.num_individuals, &|| Error::LdNoMemory {
                what: "the individuals asked for",
                values: block.num_individuals,
                bytes_per_value: size_of::<bool>(),
            })?,
        };
        for individual in individuals {
            let Some(asked_for) = asked_for_already.get_mut(*individual) else {
                return Err(Error::LdIndividualNotInTheDataset {
                    individual: *individual,
                    num_individuals: block.num_individuals,
                });
            };
            if *asked_for {
                return Err(Error::LdIndividualAskedForTwice {
                    individual: *individual,
                });
            }
            *asked_for = true;
        }
        let too_large = || Error::LdDosagesTooLarge {
            num_vars: block.num_vars,
            num_individuals,
        };
        let values = the_values_of(block.num_vars, num_individuals).ok_or_else(too_large)?;
        // The individuals each variant has a value for: the ones asked for,
        // or every individual of the block when none were named.
        let mut of_them: Vec<usize> = Vec::new();
        of_them
            .try_reserve_exact(num_individuals)
            .map_err(|_| Error::LdNoMemory {
                what: "the individuals of the dosages",
                values: num_individuals,
                bytes_per_value: size_of::<usize>(),
            })?;
        match individuals.is_empty() {
            true => of_them.extend(0..block.num_individuals),
            false => of_them.extend_from_slice(individuals),
        }
        let mut dosages = LdDosages {
            num_vars: block.num_vars,
            individuals: of_them,
            dosages: a_vector_of(
                0.0,
                values,
                &the_memory_for("the dosages", values, size_of::<f64>()),
            )?,
            called: a_vector_of(
                0.0,
                values,
                &the_memory_for("the called genotypes", values, size_of::<f64>()),
            )?,
            squares: a_vector_of(
                0.0,
                values,
                &the_memory_for("the squares of the dosages", values, size_of::<f64>()),
            )?,
            has_variance: a_vector_of(
                false,
                block.num_vars,
                &the_memory_for(
                    "the variants that have variance",
                    block.num_vars,
                    size_of::<bool>(),
                ),
            )?,
            maf: a_vector_of(
                None,
                block.num_vars,
                &the_memory_for(
                    "the major allele frequency of each variant",
                    block.num_vars,
                    size_of::<Option<f64>>(),
                ),
            )?,
        };
        if values == 0 {
            // No variant, or no individual to read at each of them: the
            // three matrices hold nothing and no genotype is read, which
            // is also what keeps the rows below from being of no allele.
            return Ok(dosages);
        }
        // The block holds its genotypes and the dosages have an individual,
        // so the block has an individual and a genotype holds an allele:
        // neither of these is `None`, and they are read as numbers that
        // cannot be 0 so that no row is cut into pieces of no allele.
        let (Some(of_a_genotype), Some(alleles_per_var)) = (
            NonZeroUsize::new(block.ploidy),
            NonZeroUsize::new(block.alleles_per_var()?),
        ) else {
            return Err(Error::GtsNotWholeGenotypes {
                num_alleles: block.gts.len(),
                ploidy: block.ploidy,
            });
        };
        // The genotypes of the individuals that were asked for, in the
        // order they were asked for, which the counts of the alleles and
        // the dosages of a variant are then read from. A call over every
        // individual reads the row of the block itself and leaves this
        // empty.
        let mut chosen = match individuals.is_empty() {
            true => Vec::new(),
            false => {
                let alleles = num_individuals
                    .checked_mul(of_a_genotype.get())
                    .ok_or_else(too_large)?;
                a_vector_of(
                    MISSING_ALLELE,
                    alleles,
                    &the_memory_for(
                        "the genotypes of the individuals asked for",
                        alleles,
                        size_of::<i8>(),
                    ),
                )?
            }
        };
        let mut counts: AlleleCounts = [0; 128];
        let rows = dosages
            .dosages
            .chunks_exact_mut(num_individuals)
            .zip(dosages.called.chunks_exact_mut(num_individuals))
            .zip(dosages.squares.chunks_exact_mut(num_individuals))
            .zip(dosages.has_variance.iter_mut())
            .zip(dosages.maf.iter_mut())
            .zip(block.gts.chunks_exact(alleles_per_var.get()));
        for (((((row, row_called), row_squares), has_variance), maf), gts) in rows {
            let genotypes = match individuals.is_empty() {
                true => gts,
                false => {
                    the_genotypes_of(gts, individuals, of_a_genotype, &mut chosen);
                    chosen.as_slice()
                }
            };
            let called_alleles = count_alleles(genotypes, &mut counts)?;
            *maf = the_major_allele_frequency(&counts, called_alleles);
            *has_variance = the_dosages_of_a_variant(
                genotypes,
                of_a_genotype,
                the_major_allele(&counts),
                row,
                row_called,
                row_squares,
            );
        }
        Ok(dosages)
    }

    /// How many variants the dosages are of.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// How many individuals each variant has a value for: the ones
    /// [`LdDosages::of_block`] was given, or every individual of the block
    /// when it was given none.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.individuals.len()
    }

    /// The variants `first..first + num_vars` of it, which the tiles of the
    /// products and the window of the filter take.
    ///
    /// # Errors
    ///
    /// When those are not variants of these dosages: the error names both
    /// numbers and how many variants there are. [`Error::LdNoMemory`] when
    /// this machine does not give the memory of the matrices of those
    /// variants, which is asked for with `try_reserve_exact` and not taken.
    pub fn rows(&self, first: usize, num_vars: usize) -> Result<LdDosages> {
        let of_other_variants = || Error::LdRowsNotInTheDosages {
            first,
            asked_for: num_vars,
            num_vars: self.num_vars,
        };
        let end = first.checked_add(num_vars).ok_or_else(of_other_variants)?;
        if end > self.num_vars {
            return Err(of_other_variants());
        }
        // The variants asked for are variants of these dosages, so both of
        // these are at most the values of one of the matrices, which is a
        // number this machine counted when they were built: neither
        // saturates.
        let from = first.saturating_mul(self.num_individuals());
        let values = num_vars.saturating_mul(self.num_individuals());
        let to = from.saturating_add(values);
        // The rows of the variants asked for, which are rows of these
        // dosages: what is not there is the error above and not a set of
        // fewer values than it says it has.
        let (Some(dosages), Some(called), Some(squares), Some(has_variance), Some(maf)) = (
            self.dosages.get(from..to),
            self.called.get(from..to),
            self.squares.get(from..to),
            self.has_variance.get(first..end),
            self.maf.get(first..end),
        ) else {
            return Err(of_other_variants());
        };
        Ok(LdDosages {
            num_vars,
            individuals: the_copy_of(
                &self.individuals,
                &the_memory_for(
                    "the individuals of the dosages",
                    self.individuals.len(),
                    size_of::<usize>(),
                ),
            )?,
            dosages: the_copy_of(
                dosages,
                &the_memory_for("the dosages", values, size_of::<f64>()),
            )?,
            called: the_copy_of(
                called,
                &the_memory_for("the called genotypes", values, size_of::<f64>()),
            )?,
            squares: the_copy_of(
                squares,
                &the_memory_for("the squares of the dosages", values, size_of::<f64>()),
            )?,
            has_variance: the_copy_of(
                has_variance,
                &the_memory_for(
                    "the variants that have variance",
                    num_vars,
                    size_of::<bool>(),
                ),
            )?,
            maf: the_copy_of(
                maf,
                &the_memory_for(
                    "the major allele frequency of each variant",
                    num_vars,
                    size_of::<Option<f64>>(),
                ),
            )?,
        })
    }

    /// Whether the called genotypes of the variant hold two dosages at
    /// least. One that does not has NaN against every variant, itself among
    /// them. It is false for a `var` that is not a variant of these
    /// dosages.
    #[must_use]
    pub fn has_variance(&self, var: usize) -> bool {
        self.has_variance.get(var).copied().unwrap_or(false)
    }

    /// The dosage of each individual at the variant, in the order the
    /// individuals were given, and `None` for a genotype with an allele
    /// missing. It is `None` for a `var` that is not a variant of these
    /// dosages.
    #[must_use]
    pub fn dosages(&self, var: usize) -> Option<impl Iterator<Item = Option<u8>> + '_> {
        if var >= self.num_vars {
            return None;
        }
        let from = var.checked_mul(self.num_individuals())?;
        let to = from.checked_add(self.num_individuals())?;
        let dosages = self.dosages.get(from..to)?;
        let called = self.called.get(from..to)?;
        Some(dosages.iter().zip(called).map(|(dosage, called)| {
            // The matrix of the called genotypes holds a 1 where the
            // genotype was called and a 0 where it was not, both written
            // here, so this asks which of the two it is.
            (*called > 0.0).then(|| the_dosage_of(*dosage))
        }))
    }

    /// The largest of the counts of the alleles over the called alleles of
    /// the variant, the major allele frequency of `docs/specs/filters.md`,
    /// over the individuals this was built with. `None` for a variant with
    /// no called allele, and for a `var` that is not a variant of these
    /// dosages.
    #[must_use]
    pub fn maf(&self, var: usize) -> Option<f64> {
        self.maf.get(var).copied().flatten()
    }
}

/// The r² of every variant of `a` against every variant of `b`, one value
/// for each pair.
///
/// `out` holds the variants of `a` one row after another, `b.num_vars()`
/// values in each, and the value of the row `i` and the column `j` is the
/// r² of the variant `i` of `a` and the variant `j` of `b`. A pair that
/// has no r² is NaN, which is what "What it gives" of `docs/specs/ld.md`
/// gives a pair whose two variants were called in no individual together,
/// and one where the dosages of either variant are all the same among the
/// individuals called at both. The r² of a variant against itself is 1
/// when it has two dosages at least and NaN when it has not.
///
/// The six sums of every pair come from six products of the three
/// matrices of the two sets, and four of them are enough when `a` and `b`
/// are the same dosages given as one reference, which is how a tile of
/// the matrix of r² against itself is asked for: the pair of the variants
/// i and j holds the same two variants as the pair of j and i, so n and
/// Σxy are the same for both and Σy and Σyy are the Σx and Σxx of the
/// pair the other way round.
///
/// # Errors
///
/// [`Error::LdR2OfAnotherSize`] when `out` does not hold one value for
/// each pair, [`Error::LdDosagesOfOtherIndividuals`] when `a` and `b`
/// were not built over the same individuals of the block in the same
/// order, and [`Error::LdLinalg`] when a product could not be worked
/// out.
pub fn r2_between(a: &LdDosages, b: &LdDosages, out: &mut [f64]) -> Result<()> {
    if let Some(problem) = the_individuals_that_differ(&a.individuals, &b.individuals) {
        return Err(Error::LdDosagesOfOtherIndividuals { problem });
    }
    let num_values = out.len();
    if a.num_vars.checked_mul(b.num_vars) != Some(num_values) {
        return Err(Error::LdR2OfAnotherSize {
            num_values,
            num_vars_of_a: a.num_vars,
            num_vars_of_b: b.num_vars,
        });
    }
    if num_values == 0 {
        // One of the two sets has no variant, so there is no pair and
        // `out` holds nothing.
        return Ok(());
    }
    if a.num_individuals() == 0 {
        // The sums of a pair run over the individuals both of its variants
        // were called in and there is no individual, so every pair has an
        // n of 0 and no r². The products are not taken, since the linear
        // algebra sums over a dimension of 1 at least. No reader of popnei
        // gives such dosages: a block of variants and no individual holds
        // no genotype, which `of_block` refuses.
        out.fill(f64::NAN);
        return Ok(());
    }
    let sums = TheSumsOfThePairs::of(a, b, num_values)?;
    match &sums.of_the_second_set {
        TheSumsOfTheSecondSet::OfTheirOwn { of_b, squares_of_b } => {
            the_r2_of_every_pair(&sums, of_b, squares_of_b, out);
        }
        TheSumsOfTheSecondSet::TheOtherWayRound => {
            the_r2_of_a_set_against_itself(&sums, a.num_vars, out);
        }
    }
    Ok(())
}

/// Writes into `out` the r² of every pair of two sets of variants, whose
/// Σy and Σyy lie where the other four sums of the pair do.
///
/// All seven buffers hold one value for each pair, in the order of the
/// pairs.
fn the_r2_of_every_pair(
    sums: &TheSumsOfThePairs,
    of_b: &[f64],
    squares_of_b: &[f64],
    out: &mut [f64],
) {
    let values = out
        .iter_mut()
        .zip(&sums.num_individuals)
        .zip(&sums.products)
        .zip(&sums.of_a)
        .zip(of_b)
        .zip(&sums.squares_of_a)
        .zip(squares_of_b);
    for ((((((r2, individuals), products), of_a), of_b), squares_of_a), squares_of_b) in values {
        *r2 = the_r2_of_a_pair(
            *individuals,
            *products,
            *of_a,
            *of_b,
            *squares_of_a,
            *squares_of_b,
        );
    }
}

/// Writes into `out` the r² of every pair of one set of `num_vars`
/// variants against itself, whose Σy and Σyy are the Σx and Σxx of the
/// pair of its two variants the other way round.
///
/// The pair of the variants i and j holds the same two variants as the
/// pair of j and i, so Σy of the row i is the column i of Σx, whose
/// entries lie one row apart, and Σyy of that row is the column i of Σxx.
/// The value read there is the one a product would have given, to the
/// bit, because the sums are whole numbers below 2^53, which
/// [`TheSumsOfTheSecondSet::TheOtherWayRound`] says why.
///
/// Every buffer holds `num_vars` rows of `num_vars` values, and `num_vars`
/// is 1 at least: a set of no variant has no pair, and the caller writes
/// nothing for it.
fn the_r2_of_a_set_against_itself(sums: &TheSumsOfThePairs, num_vars: usize, out: &mut [f64]) {
    debug_assert_eq!(
        Some(out.len()),
        num_vars.checked_mul(num_vars),
        "the r² of a set of {num_vars} variants against itself was given a buffer of {} values",
        out.len()
    );
    let rows = out
        .chunks_exact_mut(num_vars)
        .zip(sums.num_individuals.chunks_exact(num_vars))
        .zip(sums.products.chunks_exact(num_vars))
        .zip(sums.of_a.chunks_exact(num_vars))
        .zip(sums.squares_of_a.chunks_exact(num_vars))
        .enumerate();
    for (variant, ((((r2_of_the_row, individuals), products), of_a), squares_of_a)) in rows {
        let of_b = sums.of_a.iter().skip(variant).step_by(num_vars);
        let squares_of_b = sums.squares_of_a.iter().skip(variant).step_by(num_vars);
        let values = r2_of_the_row
            .iter_mut()
            .zip(individuals)
            .zip(products)
            .zip(of_a)
            .zip(of_b)
            .zip(squares_of_a)
            .zip(squares_of_b);
        for ((((((r2, individuals), products), of_a), of_b), squares_of_a), squares_of_b) in values
        {
            *r2 = the_r2_of_a_pair(
                *individuals,
                *products,
                *of_a,
                *of_b,
                *squares_of_a,
                *squares_of_b,
            );
        }
    }
}

/// How the individuals of the two sets of dosages differ, and `None` when
/// both hold the same individuals of the block in the same order.
fn the_individuals_that_differ(of_a: &[usize], of_b: &[usize]) -> Option<TheIndividualsThatDiffer> {
    if of_a.len() != of_b.len() {
        return Some(TheIndividualsThatDiffer::NotAsMany {
            of_a: of_a.len(),
            of_b: of_b.len(),
        });
    }
    of_a.iter()
        .zip(of_b)
        .enumerate()
        .find(|(_, (of_a, of_b))| of_a != of_b)
        .map(|(at, (of_a, of_b))| TheIndividualsThatDiffer::NotTheSame {
            at,
            of_a: *of_a,
            of_b: *of_b,
        })
}

/// The six sums of every pair of two sets of variants, each one a matrix
/// of the variants of the first set by those of the second, row after row.
///
/// They are the six products of "How it runs" of `docs/specs/ld.md`, and
/// each sum runs over the individuals both variants of the pair were
/// called in: the matrix of the called genotypes of one set holds a 0
/// where a genotype was not called, so every product it is in leaves that
/// individual out.
struct TheSumsOfThePairs {
    /// n, how many individuals both variants of the pair were called in.
    num_individuals: Vec<f64>,
    /// Σxy, the sum of the products of the two dosages of each of those
    /// individuals.
    products: Vec<f64>,
    /// Σx, the sum of the dosages of the variant of the first set over
    /// those individuals.
    of_a: Vec<f64>,
    /// Σxx, the sum of the squares of the dosages of the variant of the
    /// first set over them.
    squares_of_a: Vec<f64>,
    /// Where Σy and Σyy of each pair are, which is the sum of the dosages
    /// of the variant of the second set over those individuals and the sum
    /// of their squares.
    of_the_second_set: TheSumsOfTheSecondSet,
}

/// Where Σy and Σyy of each pair of two sets of variants are.
enum TheSumsOfTheSecondSet {
    /// Two products of their own, one value for each pair in the order of
    /// the pairs, which is what two sets that are not one set against
    /// itself need.
    OfTheirOwn {
        /// Σy, the sum of the dosages of the variant of the second set
        /// over the individuals both variants of the pair were called in.
        of_b: Vec<f64>,
        /// Σyy, the sum of the squares of those dosages.
        squares_of_b: Vec<f64>,
    },
    /// Σx and Σxx of the pair of its two variants the other way round,
    /// which is what one set of variants against itself has: the pair of
    /// the variants i and j holds the same two variants as the pair of j
    /// and i, so two of the six products are not taken.
    ///
    /// The two are the same number and not two numbers within a tolerance
    /// because every entry of the three matrices is a whole number and
    /// [`MAX_ALLELES_OF_A_VARIANT`] keeps every sum below the 2^53 an
    /// `f64` counts one by one. Each of them is a sum of the same
    /// products, and a routine of BLAS need not add those in the same
    /// order when it is given the two matrices in the other roles: over
    /// values that are not whole numbers the two orders differ in their
    /// last bit, which was measured on Accelerate on 23 September 2026.
    /// Anything else that reads these sums the other way round, a
    /// standardized dosage or a kinship, has to work that bound out for
    /// itself.
    TheOtherWayRound,
}

impl TheSumsOfThePairs {
    /// The six sums of every pair of `a` and `b`, which hold the same
    /// individuals in the same order and have `num_values` pairs between
    /// them, the variants of `a` times those of `b`, a number the caller
    /// has counted.
    ///
    /// The three matrices of a set of dosages hold one row for each
    /// variant and one column for each individual, and the sums of a pair
    /// run over the individuals, so the matrix of the second set is the
    /// operand of [`product`] with one row for each column of the result,
    /// which the routine reads that way and copies nothing for.
    ///
    /// # Errors
    ///
    /// [`Error::LdLinalg`] when a product could not be worked out, and
    /// [`Error::LdNoMemory`] when this machine did not give the memory of
    /// one of the sums.
    fn of(a: &LdDosages, b: &LdDosages, num_values: usize) -> Result<TheSumsOfThePairs> {
        let (rows, inner, cols) = (a.num_vars, a.num_individuals(), b.num_vars);
        let mut num_individuals = a_vector_of(
            0.0,
            num_values,
            &the_memory_for("n of the r²", num_values, size_of::<f64>()),
        )?;
        let mut products = a_vector_of(
            0.0,
            num_values,
            &the_memory_for("Σxy of the r²", num_values, size_of::<f64>()),
        )?;
        let mut of_a = a_vector_of(
            0.0,
            num_values,
            &the_memory_for("Σx of the r²", num_values, size_of::<f64>()),
        )?;
        let mut squares_of_a = a_vector_of(
            0.0,
            num_values,
            &the_memory_for("Σxx of the r²", num_values, size_of::<f64>()),
        )?;
        let sum_of = |of_the_variants: &[f64], of_the_others: &[f64], into: &mut [f64], sum| {
            // The three matrices of both sets hold one row for each
            // variant and one column for each individual, and the sums of
            // a pair run over the individuals, so the second set is the
            // operand with one row for each column of the result.
            let of_the_others = TheSecondOperand::ByTheColumnsOfTheResult {
                values: of_the_others,
                cols,
            };
            let of_the_variants = TheFirstOperand::ByTheRowsOfTheResult {
                values: of_the_variants,
                rows,
            };
            product(of_the_variants, inner, of_the_others, into).map_err(|source| Error::LdLinalg {
                operation: sum,
                source,
            })
        };
        sum_of(&a.called, &b.called, &mut num_individuals, "n")?;
        sum_of(&a.dosages, &b.dosages, &mut products, "Σxy")?;
        sum_of(&a.dosages, &b.called, &mut of_a, "Σx")?;
        sum_of(&a.squares, &b.called, &mut squares_of_a, "Σxx")?;
        let of_the_second_set = if std::ptr::eq(a, b) {
            // One set of variants against itself: the pair of the variants
            // i and j holds the two variants of the pair of j and i the
            // other way round, so Σy and Σyy are Σx and Σxx read that way
            // and two of the six products are not taken.
            TheSumsOfTheSecondSet::TheOtherWayRound
        } else {
            let mut of_b = a_vector_of(
                0.0,
                num_values,
                &the_memory_for("Σy of the r²", num_values, size_of::<f64>()),
            )?;
            let mut squares_of_b = a_vector_of(
                0.0,
                num_values,
                &the_memory_for("Σyy of the r²", num_values, size_of::<f64>()),
            )?;
            sum_of(&a.called, &b.dosages, &mut of_b, "Σy")?;
            sum_of(&a.called, &b.squares, &mut squares_of_b, "Σyy")?;
            TheSumsOfTheSecondSet::OfTheirOwn { of_b, squares_of_b }
        };
        Ok(TheSumsOfThePairs {
            num_individuals,
            products,
            of_a,
            squares_of_a,
            of_the_second_set,
        })
    }
}

/// The r² of one pair from its six sums, and NaN where it has none.
///
/// It is the formula of "What it gives" of `docs/specs/ld.md`, with n the
/// individuals both variants of the pair were called in and every sum
/// running over those individuals alone,
///
/// ```text
/// r² = (n·Σxy − Σx·Σy)² / ((n·Σxx − (Σx)²) · (n·Σyy − (Σy)²))
/// ```
///
/// A pair has no r² when n is 0 and when the dosages of one of its two
/// variants are all the same among those individuals, and each of those
/// leaves one of the two factors below the line at 0, which is the one
/// test made here. Neither factor is ever below 0: n times the sum of the
/// squares, less the square of the sum, is n² times the variance of the
/// dosages.
///
/// The six sums are whole numbers that an `f64` holds exactly, and so are
/// the four products taken of them here, n·Σxy, Σx·Σy, n·Σxx and (Σx)²,
/// and the two differences of "above the line" and of each spread, while
/// the individuals times the ploidy are at most
/// [`MAX_ALLELES_OF_A_VARIANT`], which [`LdDosages::of_block`] refuses a
/// block above. Three operations round after that: the square of what is
/// above the line, the product of the two spreads, and the division.
fn the_r2_of_a_pair(
    individuals: f64,
    products: f64,
    of_a: f64,
    of_b: f64,
    squares_of_a: f64,
    squares_of_b: f64,
) -> f64 {
    let above_the_line = individuals * products - of_a * of_b;
    let spread_of_a = individuals * squares_of_a - of_a * of_a;
    let spread_of_b = individuals * squares_of_b - of_b * of_b;
    if spread_of_a <= 0.0 || spread_of_b <= 0.0 {
        // The rule of the spec, written out. The division would give NaN
        // here without it: with the sums exact, a spread of 0 is a variant
        // whose dosages do not vary among the individuals of the pair, and
        // what is above the line is then 0 as well, so the pair comes out
        // 0/0. The test is what says so, and it is not there because the
        // arithmetic reaches a case the division gets wrong.
        return f64::NAN;
    }
    above_the_line * above_the_line / (spread_of_a * spread_of_b)
}

/// The r² of every pair of the variants a reader gives, with the
/// chromosome and the position of each of them.
///
/// It asks the reader for the genotypes, the chromosome and the position,
/// and reads it to its end. The reader is borrowed and not taken, so that
/// whoever built the chain of filters of the pass reads their counts from
/// it when this returns, as `docs/specs/filters.md` says; how many
/// variants the calculation took is [`R2Matrix::num_vars`].
///
/// `max_num_vars` is how many variants the calculation takes before it
/// refuses. The matrix holds one r² for each pair of them, so it grows
/// with the square of the variants, and a pass of more than that number is
/// an error and not a matrix this machine is asked for the memory of.
/// [`MAX_NUM_VARS_OF_THE_MATRIX`] is the number a Python or a TypeScript
/// user gets when they name none.
///
/// The r² of a pair is the one [`r2_between`] gives, and the matrix is the
/// same, to the bit, whatever the size of the blocks the reader gives and
/// however many threads the products run on: the six sums of a pair are
/// whole numbers that an `f64` holds exactly and they run over the
/// individuals, which no block and no tile cuts.
///
/// # Errors
///
/// [`Error::LdMaxNumVarsTooLarge`] when the matrix of `max_num_vars`
/// variants holds more values than this machine counts, which is looked at
/// before the pass; [`Error::LdTooManyVars`] when the pass gives more
/// variants than that, with both numbers and the memory the matrix would
/// have needed; [`Error::ReaderGaveNoVariants`] when the reader has no
/// variant; [`Error::FieldsNotInTheBlock`] when a block holds variants and
/// no genotypes or no position; [`Error::LdNoMemory`] when this machine
/// does not give the memory of the matrix, which is asked of it with
/// `try_reserve_exact` and not taken; what the dosages of a block and the
/// r² of two tiles refuse; and whatever the reader fails with, which is
/// given on as it is.
pub fn calc_r2_matrix<R: BlockReader + ?Sized>(
    reader: &mut R,
    max_num_vars: usize,
) -> Result<R2Matrix> {
    the_r2_matrix_in_tiles_of(reader, max_num_vars, THE_VARS_OF_A_TILE)
}

/// The r² of every pair of the variants of a set, with the chromosome and
/// the position of each of them.
///
/// It is what [`calc_r2_matrix`] gives for a pass over a reader.
#[derive(Debug)]
pub struct R2Matrix {
    /// How many variants the matrix is of.
    num_vars: usize,
    /// The r² of every pair, `num_vars` rows of `num_vars` values.
    r2: Vec<f64>,
    /// The number of the chromosome of each variant, in `chrom_table`.
    chroms: Vec<u32>,
    /// The names of the chromosomes, cloned from the reader of the pass,
    /// so that the numbers above are read after that reader is gone.
    chrom_table: ChromTable,
    /// The position of each variant, 1 based as in a VCF.
    poss: Vec<u64>,
}

impl R2Matrix {
    /// How many variants the matrix is of, which is how many the pass
    /// gave.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// The r² of every pair, `num_vars` x `num_vars`, row after row.
    ///
    /// The value of the row `i` and the column `j` is the r² of the
    /// variants `i` and `j` of the pass, in the order the reader gave
    /// them, and a pair that has no r² is NaN, which "What it gives" of
    /// `docs/specs/ld.md` defines. The two cells of a pair hold the same
    /// value, and the diagonal is 1 for a variant that has two dosages at
    /// least among its called genotypes and NaN for one that has not.
    #[must_use]
    pub fn r2(&self) -> &[f64] {
        &self.r2
    }

    /// The number of the chromosome of each variant, one for each row of
    /// the matrix, whose name is [`R2Matrix::chrom_table`].
    #[must_use]
    pub fn chroms(&self) -> &[u32] {
        &self.chroms
    }

    /// The names of the chromosomes of the pass, each with the number the
    /// rows of the matrix hold.
    #[must_use]
    pub fn chrom_table(&self) -> &ChromTable {
        &self.chrom_table
    }

    /// The position of each variant, 1 based as in a VCF, one for each row
    /// of the matrix. The distance of a pair is the difference of two of
    /// them, and a pair whose variants are on two chromosomes has none.
    #[must_use]
    pub fn poss(&self) -> &[u64] {
        &self.poss
    }

    /// The matrix and its three columns given away, which is how a binding
    /// crate hands them to its language.
    ///
    /// Section 1 of `docs/architecture.md` has a reader give a block away
    /// so that the binding hands the array of the genotypes over without
    /// copying it, and the matrix of the r² is the larger of the two: 200
    /// MB at [`MAX_NUM_VARS_OF_THE_MATRIX`], which a binding that read
    /// [`R2Matrix::r2`] would copy, and which in a browser is never given
    /// back. A caller in Rust that reads the matrix and keeps it uses the
    /// four accessors above instead.
    #[must_use]
    pub fn given_away(self) -> TheMatrixGivenAway {
        let R2Matrix {
            num_vars,
            r2,
            chroms,
            chrom_table,
            poss,
        } = self;
        TheMatrixGivenAway {
            num_vars,
            r2,
            chroms,
            chrom_table,
            poss,
        }
    }
}

/// What an [`R2Matrix`] holds, given away by value: the r² of every pair
/// and the three columns that say which variant each row of it is.
///
/// [`R2Matrix::given_away`] is where it comes from, and the doc comments of
/// the accessors of `R2Matrix` say what each of these is.
#[derive(Debug)]
pub struct TheMatrixGivenAway {
    /// How many variants the matrix is of.
    pub num_vars: usize,
    /// The r² of every pair, `num_vars` rows of `num_vars` values.
    pub r2: Vec<f64>,
    /// The number of the chromosome of each variant, in `chrom_table`.
    pub chroms: Vec<u32>,
    /// The names of the chromosomes, each with the number the rows hold.
    pub chrom_table: ChromTable,
    /// The position of each variant, 1 based as in a VCF.
    pub poss: Vec<u64>,
}

/// The matrix of every pair of the variants of `reader`, taken in tiles of
/// `vars_per_tile` variants.
///
/// [`calc_r2_matrix`] is this with the tile of the module, and the tests
/// are what give another: the matrix is the same, to the bit, whatever the
/// tile, because a sum of a pair runs over the individuals and the tiles
/// cut the variants.
///
/// # Errors
///
/// Those of [`calc_r2_matrix`].
fn the_r2_matrix_in_tiles_of<R: BlockReader + ?Sized>(
    reader: &mut R,
    max_num_vars: usize,
    vars_per_tile: usize,
) -> Result<R2Matrix> {
    // The matrix holds the square of the variants of the pass and the pass
    // takes `max_num_vars` of them at most, so this is what says that every
    // count of values below is a number this machine counts: the largest
    // number whose square a `usize` holds, 65535 in WebAssembly, where a
    // `usize` is 32 bits, and 4294967295 natively.
    if max_num_vars.checked_mul(max_num_vars).is_none() {
        return Err(Error::LdMaxNumVarsTooLarge { max_num_vars });
    }
    // A tile of no variant would take no variant of a block and the pass
    // would stand still.
    let vars_per_tile = vars_per_tile.max(1);
    let pass = the_dosages_of_the_pass(reader, max_num_vars, vars_per_tile)?;
    if pass.num_vars == 0 {
        let filters = reader.filtering_stats();
        return Err(Error::PassGaveNoVariant {
            // The filter nearest the source was given what the source
            // gave; with no filter the pass gave what the source gave,
            // which is nothing.
            num_vars_of_the_source: filters.last().map_or(0, |(_, stats)| stats.vars_processed),
            filters,
        });
    }
    let r2 = the_r2_of_the_tiles(&pass.tiles, pass.num_vars, max_num_vars)?;
    Ok(R2Matrix {
        num_vars: pass.num_vars,
        r2,
        chroms: pass.chroms,
        chrom_table: reader.chroms().clone(),
        poss: pass.poss,
    })
}

/// The dosages of the variants of a pass, in tiles, with the chromosome
/// and the position of each variant.
struct ThePassOfTheMatrix {
    /// How many variants the pass gave, which is the variants of its tiles
    /// together.
    num_vars: usize,
    /// The dosages of those variants: the tiles of the products, in the
    /// order of the variants.
    tiles: Vec<LdDosages>,
    /// The number of the chromosome of each variant, in the table of the
    /// reader.
    chroms: Vec<u32>,
    /// The position of each variant, 1 based as in a VCF.
    poss: Vec<u64>,
}

/// Reads `reader` to its end and gives the dosages of its variants in
/// tiles of `vars_per_tile` variants, with the chromosome and the position
/// of each.
///
/// # Errors
///
/// [`Error::LdTooManyVars`] as soon as the variants pass `max_num_vars`,
/// so that a source of a million variants is not read to its end to be
/// refused; [`Error::FieldsNotInTheBlock`] when a block holds variants and
/// no genotypes or no position; [`Error::BlocksDoNotFitTogether`] when a
/// block holds other individuals or another ploidy than the reader says
/// its source has; [`Error::LdNoMemory`] when this machine does not give
/// the memory of the chromosomes and the positions; what
/// [`Block::check`] and the dosages of a tile refuse; and whatever the
/// reader fails with.
fn the_dosages_of_the_pass<R: BlockReader + ?Sized>(
    reader: &mut R,
    max_num_vars: usize,
    vars_per_tile: usize,
) -> Result<ThePassOfTheMatrix> {
    // The genotypes, the chromosome and the position are what this reads,
    // so a reader over a file leaves the other columns of a variant
    // unparsed.
    reader.set_needs(Needs::GTS | Needs::CHROM_POS);
    let mut tiles =
        TheTilesOfThePass::of(vars_per_tile, reader.individuals().len(), reader.ploidy());
    let mut num_vars = 0_usize;
    let mut chroms: Vec<u32> = Vec::new();
    let mut poss: Vec<u64> = Vec::new();
    while let Some(block) = reader.next_block()? {
        block.check()?;
        let missing = (Needs::GTS | Needs::CHROM_POS).difference(block.fields());
        if !missing.is_empty() {
            return Err(Error::FieldsNotInTheBlock { fields: missing });
        }
        if block.num_individuals != tiles.num_individuals || block.ploidy != tiles.ploidy {
            return Err(Error::BlocksDoNotFitTogether {
                num_individuals: tiles.num_individuals,
                ploidy: tiles.ploidy,
                found_num_individuals: block.num_individuals,
                found_ploidy: block.ploidy,
            });
        }
        // The count is refused on the next line as soon as it passes the
        // cap, so what saturates here is a pass that no machine gave: the
        // message then names the largest number a `usize` holds.
        let with_the_block = num_vars.saturating_add(block.num_vars);
        if with_the_block > max_num_vars {
            return Err(the_variants_pass_the_cap(with_the_block, max_num_vars));
        }
        num_vars = with_the_block;
        // The block holds the two columns, which the fields above say.
        let (Some(of_its_variants), Some(at_which_they_are)) = (&block.chrom, &block.pos) else {
            return Err(Error::FieldsNotInTheBlock {
                fields: Needs::CHROM_POS,
            });
        };
        the_values_of_the_column(
            &mut chroms,
            of_its_variants,
            "the chromosome of each variant",
        )?;
        the_values_of_the_column(&mut poss, at_which_they_are, "the position of each variant")?;
        tiles.take_the_block(&block)?;
        // The block is given back before the reader is asked for the next
        // one, so the memory of two blocks is never held at once.
        drop(block);
    }
    Ok(ThePassOfTheMatrix {
        num_vars,
        tiles: tiles.done()?,
        chroms,
        poss,
    })
}

/// The tiles of a pass, which the blocks of the reader are poured into.
///
/// A tile is the dosages of `vars_per_tile` variants of the pass, the last
/// one of what is left, and the tiles are cut at multiples of that number
/// counted from the first variant of the pass and not where the blocks
/// end. So the products are taken over the same variants together whatever
/// the reader gives at a time.
///
/// Each tile is built from the genotypes of its own variants, gathered
/// from the blocks into one buffer that the whole pass shares, and holds
/// the three matrices of those variants. A pair of tiles is then a pair of
/// operands of the products as it stands, and nothing of a tile is copied
/// to take one.
struct TheTilesOfThePass {
    /// How many variants a tile holds, 1 at least.
    vars_per_tile: usize,
    /// How many individuals the reader says its source has, which is what
    /// the dosages of every tile are built over.
    num_individuals: usize,
    /// How many alleles the genotype of one individual holds.
    ploidy: usize,
    /// The genotypes of the tile being filled, variant after variant. It
    /// is taken back from each tile with the memory it has, so a pass over
    /// a million variants allocates it once.
    gts: Vec<i8>,
    /// How many variants of the tile being filled are in `gts`.
    vars_of_the_tile: usize,
    /// The tiles that are full, in the order of the variants.
    tiles: Vec<LdDosages>,
}

impl TheTilesOfThePass {
    /// The tiles of a pass over a source of `num_individuals` individuals
    /// of the ploidy `ploidy`, each of `vars_per_tile` variants.
    fn of(vars_per_tile: usize, num_individuals: usize, ploidy: usize) -> TheTilesOfThePass {
        TheTilesOfThePass {
            vars_per_tile,
            num_individuals,
            ploidy,
            gts: Vec::new(),
            vars_of_the_tile: 0,
            tiles: Vec::new(),
        }
    }

    /// Puts the variants of `block` into the tiles, building each tile as
    /// soon as its variants are all there.
    ///
    /// # Errors
    ///
    /// [`Error::LdNoMemory`] when this machine does not give the memory of
    /// the genotypes of a tile, [`Error::BlockArrayOfAnotherSize`] when
    /// the genotypes of the block are not its variants times the alleles
    /// of one variant, and what [`LdDosages::of_block`] refuses.
    fn take_the_block(&mut self, block: &Block) -> Result<()> {
        let alleles_per_var = block.alleles_per_var()?;
        // How many variants of the block are in a tile already.
        let mut given = 0_usize;
        while given < block.num_vars {
            // The tile has room for a variant and the block has one left,
            // so both of these are 1 at least and the loop moves on.
            let room = self.vars_per_tile.saturating_sub(self.vars_of_the_tile);
            let taken = room.min(block.num_vars.saturating_sub(given));
            // `Block::check` has passed, so the genotypes of the block are
            // its variants times the alleles of one variant, a number this
            // machine counted: neither of these saturates.
            let from = given.saturating_mul(alleles_per_var);
            let to = from.saturating_add(taken.saturating_mul(alleles_per_var));
            let Some(genotypes) = block.gts.get(from..to) else {
                return Err(Error::BlockArrayOfAnotherSize {
                    array: "gts",
                    found: block.gts.len(),
                    expected: to,
                });
            };
            self.gts
                .try_reserve(genotypes.len())
                .map_err(|_| Error::LdNoMemory {
                    what: "the genotypes of a tile",
                    values: genotypes.len(),
                    bytes_per_value: size_of::<i8>(),
                })?;
            self.gts.extend_from_slice(genotypes);
            self.vars_of_the_tile = self.vars_of_the_tile.saturating_add(taken);
            given = given.saturating_add(taken);
            if self.vars_of_the_tile == self.vars_per_tile {
                self.the_tile_is_full()?;
            }
        }
        Ok(())
    }

    /// The tiles of the pass, with the last one, which holds the variants
    /// left over from the tile before it.
    ///
    /// # Errors
    ///
    /// What [`LdDosages::of_block`] refuses.
    fn done(mut self) -> Result<Vec<LdDosages>> {
        if self.vars_of_the_tile > 0 {
            self.the_tile_is_full()?;
        }
        Ok(self.tiles)
    }

    /// Builds the dosages of the variants in the buffer and keeps them as
    /// the next tile, with the buffer left empty for the tile after it.
    ///
    /// # Errors
    ///
    /// What [`LdDosages::of_block`] refuses, and [`Error::LdNoMemory`]
    /// when this machine does not give the memory of the tiles.
    fn the_tile_is_full(&mut self) -> Result<()> {
        let mut block = Block {
            num_vars: self.vars_of_the_tile,
            num_individuals: self.num_individuals,
            ploidy: self.ploidy,
            gts: std::mem::take(&mut self.gts),
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        };
        let tile = LdDosages::of_block(&block, &[])?;
        // The buffer is taken back with the memory it has and the tile
        // after this one is gathered into it.
        self.gts = std::mem::take(&mut block.gts);
        self.gts.clear();
        self.vars_of_the_tile = 0;
        self.tiles.try_reserve(1).map_err(|_| Error::LdNoMemory {
            what: "the tiles of the products",
            values: self.tiles.len(),
            bytes_per_value: size_of::<LdDosages>(),
        })?;
        self.tiles.push(tile);
        Ok(())
    }
}

/// Adds the values of a column of a block to the ones the blocks before it
/// gave.
///
/// # Errors
///
/// [`Error::LdNoMemory`], which `what` names the column of, when this
/// machine does not give the memory of the values.
fn the_values_of_the_column<T: Copy>(
    of_the_pass: &mut Vec<T>,
    of_the_block: &[T],
    what: &'static str,
) -> Result<()> {
    of_the_pass
        .try_reserve(of_the_block.len())
        .map_err(|_| Error::LdNoMemory {
            what,
            values: of_the_block.len(),
            bytes_per_value: size_of::<T>(),
        })?;
    of_the_pass.extend_from_slice(of_the_block);
    Ok(())
}

/// The error of a pass of `num_vars` variants where `max_num_vars` were
/// allowed, with the memory the matrix of those variants would have
/// needed.
fn the_variants_pass_the_cap(num_vars: usize, max_num_vars: usize) -> Error {
    // The bytes of a matrix that was never asked for: they are counted in
    // a `u64` and not in a `usize`, so that the number in the message is
    // the right one in WebAssembly too, where a `usize` is 32 bits and the
    // matrix of 23171 variants is already more bytes than one counts.
    let of_a_variant = u64::try_from(num_vars).unwrap_or(u64::MAX);
    let bytes = of_a_variant.saturating_mul(of_a_variant).saturating_mul(8);
    Error::LdTooManyVars {
        num_vars,
        max_num_vars,
        bytes,
    }
}

/// The r² of every pair of the variants of the tiles, `num_vars` rows of
/// `num_vars` values, row after row.
///
/// The pairs are taken tile pair by tile pair, and only the pairs of tiles
/// from the diagonal up: r² is the same whichever variant of a pair comes
/// first, so the products of a pair of tiles are taken once and their
/// values are written into the two halves of the matrix. A tile against
/// itself is given to [`r2_between`] as one reference twice, which is what
/// makes it take the four products of a set against itself and not six.
///
/// `max_num_vars` is what the pass was allowed, which the caller has found
/// to have a square this machine counts.
///
/// # Errors
///
/// [`Error::LdMaxNumVarsTooLarge`] when the square of the variants of the
/// pass is not a number this machine counts, which the caller's check makes
/// unreachable; [`Error::LdNoMemory`] when this machine does not give the
/// memory of the matrix or of the r² of one pair of tiles; and what
/// [`r2_between`] refuses.
fn the_r2_of_the_tiles(
    tiles: &[LdDosages],
    num_vars: usize,
    max_num_vars: usize,
) -> Result<Vec<f64>> {
    // The caller has refused a `max_num_vars` whose square is not a number
    // this machine counts and the pass gave at most that many variants, so
    // the square below is there. The error names the number the user wrote
    // and not the variants of the pass, which is the number they would
    // lower: a pass of more variants than `max_num_vars` was stopped before
    // this.
    let values = num_vars
        .checked_mul(num_vars)
        .ok_or(Error::LdMaxNumVarsTooLarge { max_num_vars })?;
    let mut matrix = a_vector_of(
        f64::NAN,
        values,
        &the_memory_for(
            "the matrix of the r² of every pair",
            values,
            size_of::<f64>(),
        ),
    )?;
    // The r² of one pair of tiles, which every pair of them is written
    // into: the largest tile against itself, 512 KB at the 256 variants of
    // `THE_VARS_OF_A_TILE`.
    let of_the_largest = tiles.iter().map(LdDosages::num_vars).max().unwrap_or(0);
    let values = of_the_largest.saturating_mul(of_the_largest);
    let mut of_the_pair = a_vector_of(
        0.0,
        values,
        &the_memory_for("the r² of a pair of tiles", values, size_of::<f64>()),
    )?;
    let mut first_row = 0_usize;
    for (of_a, tile_a) in tiles.iter().enumerate() {
        let mut first_col = first_row;
        for tile_b in tiles.iter().skip(of_a) {
            let (rows, cols) = (tile_a.num_vars(), tile_b.num_vars());
            let values = rows.saturating_mul(cols);
            let Some(of_the_pair) = of_the_pair.get_mut(..values) else {
                return Err(Error::LdR2OfAnotherSize {
                    num_values: values,
                    num_vars_of_a: rows,
                    num_vars_of_b: cols,
                });
            };
            // On the diagonal `tile_a` and `tile_b` are the same tile, and
            // the two arguments are then one reference given twice, which
            // is what the four products of a set against itself are taken
            // on.
            r2_between(tile_a, tile_b, of_the_pair)?;
            write_the_pair_of_tiles(
                of_the_pair,
                (first_row, first_col),
                (rows, cols),
                &mut matrix,
                num_vars,
            )?;
            first_col = first_col.saturating_add(cols);
        }
        first_row = first_row.saturating_add(tile_a.num_vars());
    }
    Ok(matrix)
}

/// Writes the r² of a pair of tiles into the matrix of the pass, at the
/// rows of the first tile and the columns of the second and at the cells
/// the other way round.
///
/// `of_the_pair` holds `rows` rows, one for each variant of the first
/// tile, and `cols` values in each, the r² of that variant against the
/// variants of the second tile; the matrix holds `num_vars` rows of
/// `num_vars` values, and the two tiles begin at the variants
/// `first_row` and `first_col` of the pass. A pair of tiles on the
/// diagonal writes the square of its own variants and nothing more: its
/// two halves are already in `of_the_pair`.
///
/// # Errors
///
/// [`Error::LdRowsNotInTheDosages`] when the cells of the pair are not
/// cells of the matrix, which is a defect of the tiling and not anything a
/// caller of the crate wrote.
fn write_the_pair_of_tiles(
    of_the_pair: &[f64],
    (first_row, first_col): (usize, usize),
    (rows, cols): (usize, usize),
    matrix: &mut [f64],
    num_vars: usize,
) -> Result<()> {
    let of_other_cells = |first: usize, asked_for: usize| Error::LdRowsNotInTheDosages {
        first,
        asked_for,
        num_vars,
    };
    if cols == 0 {
        return Ok(());
    }
    // Every row of the tile is a run of the matrix: the values of one
    // variant of the first tile against the variants of the second lie
    // side by side in the row of that variant.
    for (row, values) in of_the_pair.chunks_exact(cols).enumerate() {
        // The tile is inside the matrix, whose values are the square of
        // the variants of the pass, a number this machine counted: none of
        // these saturates.
        let from = first_row
            .saturating_add(row)
            .saturating_mul(num_vars)
            .saturating_add(first_col);
        let Some(into) = matrix.get_mut(from..from.saturating_add(cols)) else {
            return Err(of_other_cells(first_col, cols));
        };
        into.copy_from_slice(values);
    }
    if first_row == first_col {
        // The tile is on the diagonal, so it is square and holds the two
        // halves of its own variants, which the rows above wrote.
        return Ok(());
    }
    // And the same values the other way round: the cells of a variant of
    // the second tile lie side by side in the row of that variant, where
    // they are the column of `of_the_pair` of that variant.
    for column in 0..cols {
        let from = first_col
            .saturating_add(column)
            .saturating_mul(num_vars)
            .saturating_add(first_row);
        let Some(cells) = matrix.get_mut(from..from.saturating_add(rows)) else {
            return Err(of_other_cells(first_row, rows));
        };
        for (cell, values) in cells.iter_mut().zip(of_the_pair.chunks_exact(cols)) {
            let Some(value) = values.get(column) else {
                return Err(of_other_cells(first_row, rows));
            };
            *cell = *value;
        }
    }
    Ok(())
}

/// `values` copies of `value`, or the error that `not_given` builds when
/// this machine did not give the memory for them.
///
/// The memory is asked for with `try_reserve_exact`, which gives it back
/// as an error where `vec![value; values]` would end the process. It is
/// also what refuses a number of values whose bytes this machine does not
/// count: a `usize` is 32 bits in WebAssembly, where 2^29 values of 8
/// bytes are already more than one holds.
fn a_vector_of<T: Clone>(
    value: T,
    values: usize,
    not_given: &impl Fn() -> Error,
) -> Result<Vec<T>> {
    let mut vector: Vec<T> = Vec::new();
    vector.try_reserve_exact(values).map_err(|_| not_given())?;
    vector.resize(values, value);
    Ok(vector)
}

/// The error of `values` values of `bytes_per_value` bytes that this machine
/// did not give the memory for, which `what` names.
///
/// `bytes_per_value` is `size_of` of the value the vector holds, and it is
/// given at the call and not worked out here because the type of the vector
/// is the one its value is inferred from: the r² is held in `f64`, the
/// genotypes of the window of a population in `i8`, and a variant of that
/// window in a chromosome and a position.
fn the_memory_for(what: &'static str, values: usize, bytes_per_value: usize) -> impl Fn() -> Error {
    move || Error::LdNoMemory {
        what,
        values,
        bytes_per_value,
    }
}

/// How many values a matrix of `num_vars` variants of `num_individuals`
/// individuals holds, or `None` when they are more than the linear algebra
/// counts in.
fn the_values_of(num_vars: usize, num_individuals: usize) -> Option<usize> {
    num_vars
        .checked_mul(num_individuals)
        .filter(|values| *values <= MAX_VALUES_OF_THE_DOSAGES)
}

/// A copy of `values`, with its memory asked of this machine with
/// `try_reserve_exact` and not taken, as `docs/specs/linalg.md` asks for
/// the workspace of the eigendecomposition.
///
/// # Errors
///
/// What `not_given` gives, when the machine does not give the memory.
fn the_copy_of<T: Copy>(values: &[T], not_given: &impl Fn() -> Error) -> Result<Vec<T>> {
    let mut copy: Vec<T> = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| not_given())?;
    copy.extend_from_slice(values);
    Ok(copy)
}

/// The genotypes of the individuals of `individuals`, in the order they are
/// given, written into `chosen`.
///
/// `gts` is one row of the genotypes of a block, the alleles of one
/// individual after those of the individual before it, and `chosen` holds
/// one genotype for each index. [`LdDosages::of_block`] refuses an index
/// that is not an individual of the block before it reads a variant, so
/// every genotype is found; one that was not is written as a genotype with
/// no allele called, since `chosen` is kept from one variant to the next
/// and would otherwise hold the genotype of the variant before it.
fn the_genotypes_of(
    gts: &[i8],
    individuals: &[usize],
    of_a_genotype: NonZeroUsize,
    chosen: &mut [i8],
) {
    for (genotype, individual) in chosen
        .chunks_exact_mut(of_a_genotype.get())
        .zip(individuals)
    {
        // An individual of the row is below the individuals of the block,
        // and the alleles of the block's own row are those individuals
        // times the ploidy, a number the block was built with: neither of
        // these saturates.
        let from = individual.saturating_mul(of_a_genotype.get());
        let to = from.saturating_add(of_a_genotype.get());
        match gts.get(from..to) {
            Some(theirs) => {
                for (allele, theirs) in genotype.iter_mut().zip(theirs) {
                    *allele = *theirs;
                }
            }
            None => genotype.fill(MISSING_ALLELE),
        }
    }
}

/// The three values of each individual at one variant, and whether the
/// called genotypes of it hold two dosages at least.
///
/// `genotypes` holds one genotype of `of_a_genotype` alleles for each
/// individual, and `row`, `called` and `squares` one value each: the dosage
/// of the genotype, how many of its alleles are not `major`, with a 0 for a
/// genotype that has an allele missing; a 1 where the genotype was called
/// and a 0 where it was not; and the square of the dosage.
///
/// The ploidies 1 to 4 each get the loop with the length of a genotype
/// written into it, because with that length a constant the compiler reads
/// the genotypes of several individuals at once, and with it a number the
/// dataset carries it reads one allele at a time. Every other ploidy takes
/// [`the_dosages_of_any_ploidy`], which is the same body with the length
/// read from the dataset. The arms give the same values: the dosage is a
/// count of alleles and the missing flag is a boolean.
fn the_dosages_of_a_variant(
    genotypes: &[i8],
    of_a_genotype: NonZeroUsize,
    major: i8,
    row: &mut [f64],
    called: &mut [f64],
    squares: &mut [f64],
) -> bool {
    match of_a_genotype.get() {
        1 => the_dosages_of_a_ploidy_of::<1>(genotypes, major, row, called, squares),
        2 => the_dosages_of_a_ploidy_of::<2>(genotypes, major, row, called, squares),
        3 => the_dosages_of_a_ploidy_of::<3>(genotypes, major, row, called, squares),
        4 => the_dosages_of_a_ploidy_of::<4>(genotypes, major, row, called, squares),
        of_a_genotype => {
            the_dosages_of_any_ploidy(genotypes, of_a_genotype, major, row, called, squares)
        }
    }
}

/// The three values of each individual at one variant whose genotypes hold
/// `OF_A_GENOTYPE` alleles, which is the body of
/// [`the_dosages_of_a_variant`] with the length of a genotype known when
/// the code is compiled.
///
/// Neither the genotype with an allele missing nor the two dosages of the
/// variant branch: the first is a choice between two values the compiler
/// takes lane by lane, and the second is the lowest and the highest dosage
/// of the called genotypes, which differ exactly when two called genotypes
/// differ. A genotype that was not called leaves both alone, since it
/// offers [`u8::MAX`] to the lowest and 0 to the highest, and a variant of
/// no called genotype ends with a lowest of [`u8::MAX`] and a highest of 0,
/// which is why the answer is the strict comparison and not an inequality.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the dosage counts the alleles of one genotype, which are the ploidy, and \
              `LdDosages::of_block` refuses a ploidy above the 255 a u8 holds"
)]
fn the_dosages_of_a_ploidy_of<const OF_A_GENOTYPE: usize>(
    genotypes: &[i8],
    major: i8,
    row: &mut [f64],
    called: &mut [f64],
    squares: &mut [f64],
) -> bool {
    let (genotypes, _) = genotypes.as_chunks::<OF_A_GENOTYPE>();
    let mut lowest = u8::MAX;
    let mut highest = 0_u8;
    let values = row
        .iter_mut()
        .zip(called.iter_mut())
        .zip(squares.iter_mut())
        .zip(genotypes);
    for (((dosage_of, called_of), square_of), genotype) in values {
        let mut dosage = 0_u8;
        let mut missing = 0_u8;
        for allele in genotype {
            dosage += u8::from(*allele != major);
            missing |= u8::from(*allele == MISSING_ALLELE);
        }
        // A genotype with an allele missing has no dosage: it is a 0 in the
        // dosages and in their squares, and a 0 in the called genotypes
        // takes it out of every sum of every pair its variant is in.
        let value = if missing == 0 { f64::from(dosage) } else { 0.0 };
        *dosage_of = value;
        *called_of = if missing == 0 { 1.0 } else { 0.0 };
        *square_of = value * value;
        lowest = lowest.min(if missing == 0 { dosage } else { u8::MAX });
        highest = highest.max(if missing == 0 { dosage } else { 0 });
    }
    highest > lowest
}

/// The three values of each individual at one variant whose genotypes hold
/// `of_a_genotype` alleles, a length the dataset carries.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the dosage counts the alleles of one genotype, which are the ploidy, and \
              `LdDosages::of_block` refuses a ploidy above the 255 a u8 holds"
)]
fn the_dosages_of_any_ploidy(
    genotypes: &[i8],
    of_a_genotype: usize,
    major: i8,
    row: &mut [f64],
    called: &mut [f64],
    squares: &mut [f64],
) -> bool {
    let mut lowest = u8::MAX;
    let mut highest = 0_u8;
    let values = row
        .iter_mut()
        .zip(called.iter_mut())
        .zip(squares.iter_mut())
        .zip(genotypes.chunks_exact(of_a_genotype));
    for (((dosage_of, called_of), square_of), genotype) in values {
        let mut dosage = 0_u8;
        let mut missing = 0_u8;
        for allele in genotype {
            dosage += u8::from(*allele != major);
            missing |= u8::from(*allele == MISSING_ALLELE);
        }
        let value = if missing == 0 { f64::from(dosage) } else { 0.0 };
        *dosage_of = value;
        *called_of = if missing == 0 { 1.0 } else { 0.0 };
        *square_of = value * value;
        lowest = lowest.min(if missing == 0 { dosage } else { u8::MAX });
        highest = highest.max(if missing == 0 { dosage } else { 0 });
    }
    highest > lowest
}

/// The dosage that an entry of the matrix of the dosages holds.
///
/// # Panics
///
/// In a build with the debug assertions on, when the value is not a whole
/// number from 0 to the 255 a dosage is held in, which is what
/// [`the_dosages_of_a_variant`] writes into that matrix and nothing else
/// does.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "an entry of the matrix is a whole number from 0 to the ploidy, written by \
              `the_dosages_of_a_variant`, and `LdDosages::of_block` refuses a ploidy above \
              the 255 a u8 holds"
)]
fn the_dosage_of(value: f64) -> u8 {
    debug_assert!(
        value.is_finite() && (0.0..=255.0).contains(&value) && value.fract() == 0.0,
        "the matrix of the dosages holds {value}, which is not a dosage"
    );
    value as u8
}

/// What the benchmark `r2_matrix` calls to take the matrix with a tile of
/// its own choosing.
///
/// How many variants a tile holds is [`THE_VARS_OF_A_TILE`], which is
/// private and which [`calc_r2_matrix`] passes on: nothing outside this
/// module names another, and a benchmark is a crate of its own. The
/// performance review of `docs/plans/ld.md` is what settles that number,
/// and it settles it by timing the matrix at several tiles over one file,
/// so this module is behind the cargo feature `bench-internals`, which is
/// off by default and which nothing of popnei's own builds turn on, as
/// `variant::bench_internals` is.
///
/// It holds one wrapper, over the private function [`calc_r2_matrix`]
/// itself calls, so that what the benchmark times is the code the library
/// runs.
#[cfg(feature = "bench-internals")]
#[doc(hidden)]
pub mod bench_internals {
    use super::{R2Matrix, THE_VARS_OF_A_TILE, the_r2_matrix_in_tiles_of as r2_matrix_in_tiles_of};

    /// How many variants a tile holds when the library chooses, which is
    /// what [`calc_r2_matrix`](super::calc_r2_matrix) passes on. The
    /// benchmark takes it as the tile it runs at when its command line
    /// names none, so the number it reports by default is the number the
    /// library uses and cannot drift from it.
    pub const VARS_OF_A_TILE: usize = THE_VARS_OF_A_TILE;
    use crate::block::BlockReader;
    use crate::error::Result;

    /// The r² of every pair of the variants `reader` gives, with the
    /// products taken in tiles of `vars_per_tile` variants, which is
    /// `the_r2_matrix_in_tiles_of` of this module.
    ///
    /// [`calc_r2_matrix`](super::calc_r2_matrix) is this with the tile the
    /// library chooses, and the matrix is the same whatever the tile: the
    /// tiles cut the variants and every sum of a pair runs over the
    /// individuals.
    ///
    /// # Errors
    ///
    /// Those of [`calc_r2_matrix`](super::calc_r2_matrix).
    pub fn the_r2_matrix_in_tiles_of<R: BlockReader + ?Sized>(
        reader: &mut R,
        max_num_vars: usize,
        vars_per_tile: usize,
    ) -> Result<R2Matrix> {
        r2_matrix_in_tiles_of(reader, max_num_vars, vars_per_tile)
    }
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;
    use std::path::{Path, PathBuf};

    use super::{
        LdDosages, MAX_ALLELES_OF_A_VARIANT, MAX_NUM_VARS_OF_THE_MATRIX, MAX_PLOIDY_OF_THE_DOSAGES,
        MAX_VALUES_OF_THE_DOSAGES, R2Matrix, THE_VARS_OF_A_TILE, TheIndividualsThatDiffer,
        TheSumsOfThePairs, TheSumsOfTheSecondSet, a_vector_of, calc_r2_matrix, r2_between,
        the_dosages_of_a_variant, the_dosages_of_any_ploidy, the_genotypes_of, the_memory_for,
        the_r2_matrix_in_tiles_of, the_values_of,
    };
    use popnei_linalg::Error as LinalgError;

    use crate::block::{Block, BlockReader};
    use crate::error::Error;
    use crate::filters::FilteringStats;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{ChromTable, MISSING_ALLELE, Needs};

    /// The allele that was not called, which shortens the tables of
    /// genotypes below.
    const M: i8 = MISSING_ALLELE;

    /// The five variants of six diploid individuals of the worked example
    /// of "How it is verified" of `docs/specs/ld.md`, each row the alleles
    /// of one individual after those of the individual before it.
    ///
    /// Their dosages are the table of that example: the major allele of v1
    /// and of v5 is allele 0, which the tie between the two alleles gives
    /// to the lower numbered one, of v2 it is allele 1, and v4 has one
    /// genotype in every individual and so no variance.
    const THE_WORKED_EXAMPLE: [[i8; 12]; 5] = [
        // v1 0/0 0/0 0/1 0/1 1/1 1/1
        [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1],
        // v2 0/0 0/1 0/1 1/1 1/1 1/1
        [0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1],
        // v3 0/0 0/0 0/0 0/1 ./. 1/1
        [0, 0, 0, 0, 0, 0, 0, 1, M, M, 1, 1],
        // v4 0/0 0/0 0/0 0/0 0/0 0/0
        [0; 12],
        // v5 0/1 1/1 0/0 0/1 1/1 0/0
        [0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0],
    ];

    /// A block of the variants given, of `num_individuals` individuals of
    /// the ploidy `ploidy`, with the genotypes and no column: the fields
    /// the dosages ask their reader for.
    fn block_of(variants: &[&[i8]], num_individuals: usize, ploidy: usize) -> Block {
        let mut gts = Vec::new();
        for row in variants {
            gts.extend_from_slice(row);
        }
        Block {
            num_vars: variants.len(),
            num_individuals,
            ploidy,
            gts,
            chrom: None,
            pos: None,
            id: None,
            alleles: None,
            qual: None,
        }
    }

    /// The block of the worked example.
    fn the_worked_example() -> Block {
        let rows: Vec<&[i8]> = THE_WORKED_EXAMPLE
            .iter()
            .map(|row| row.as_slice())
            .collect();
        block_of(&rows, 6, 2)
    }

    /// The dosages of the variant, one value for each individual.
    fn dosages_of(dosages: &LdDosages, var: usize) -> Vec<Option<u8>> {
        dosages
            .dosages(var)
            .unwrap_or_else(|| panic!("the variant {var} has no dosages"))
            .collect()
    }

    /// That two rows of one of the matrices hold the same values.
    ///
    /// The entries of the three matrices are whole numbers below 2^53,
    /// which an `f64` holds exactly, so they are compared as they are and
    /// a value that is not the one expected is not a rounding.
    #[expect(
        clippy::float_cmp,
        reason = "the entries of the three matrices are whole numbers below 2^53, which an f64 holds exactly"
    )]
    fn assert_the_values_are(found: &[f64], expected: &[f64], what: &str) {
        assert_eq!(
            found.len(),
            expected.len(),
            "{what}: the values are not as many"
        );
        for (at, (found, expected)) in found.iter().zip(expected).enumerate() {
            assert!(
                *found == *expected,
                "{what}: the value {at} is {found} and not {expected}"
            );
        }
    }

    /// That a major allele frequency is the one expected, which is the
    /// division of two counts that the code under test makes, so the two
    /// sides round the same way and are compared as they are.
    #[expect(
        clippy::float_cmp,
        reason = "both sides are one division of the same two counts, which rounds the same way on both"
    )]
    fn assert_the_number_is(found: f64, expected: f64, what: &str) {
        assert!(
            found == expected,
            "{what}: it is {found} and not {expected}"
        );
    }

    #[test]
    fn the_dosages_of_the_worked_example_are_the_ones_of_the_spec() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        assert_eq!((dosages.num_vars(), dosages.num_individuals()), (5, 6));
        let of_the_spec: [[Option<u8>; 6]; 5] = [
            [Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)],
            [Some(2), Some(1), Some(1), Some(0), Some(0), Some(0)],
            [Some(0), Some(0), Some(0), Some(1), None, Some(2)],
            [Some(0), Some(0), Some(0), Some(0), Some(0), Some(0)],
            [Some(1), Some(2), Some(0), Some(1), Some(2), Some(0)],
        ];
        for (var, expected) in of_the_spec.iter().enumerate() {
            assert_eq!(
                dosages_of(&dosages, var),
                expected.to_vec(),
                "the variant {var}"
            );
        }
        assert!(dosages.dosages(5).is_none(), "a sixth variant has dosages");
    }

    #[test]
    fn the_three_matrices_are_the_dosages_the_called_genotypes_and_the_squares() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        // v3, 0/0 0/0 0/0 0/1 ./. 1/1, whose fifth individual is the one
        // genotype of the example that was not called.
        let third = 12..18;
        assert_the_values_are(
            dosages
                .dosages
                .get(third.clone())
                .expect("the dosages of v3"),
            &[0.0, 0.0, 0.0, 1.0, 0.0, 2.0],
            "the dosages of v3",
        );
        assert_the_values_are(
            dosages.called.get(third.clone()).expect("the called of v3"),
            &[1.0, 1.0, 1.0, 1.0, 0.0, 1.0],
            "the called genotypes of v3",
        );
        assert_the_values_are(
            dosages.squares.get(third).expect("the squares of v3"),
            &[0.0, 0.0, 0.0, 1.0, 0.0, 4.0],
            "the squares of the dosages of v3",
        );
    }

    /// The alleles repeat with a period of 11, which none of the four
    /// ploidies divides, so a missing allele falls in every position of a
    /// genotype, and the run of four major alleles the period opens with
    /// gives a genotype of the dosage 0 at each of them. The four
    /// assertions at the end check that the fixture reaches the three cases
    /// the two loops could disagree on: a genotype with an allele missing,
    /// a genotype of the major allele only, a variant of two dosages at
    /// least, and a variant of one called dosage, which is the last one.
    #[test]
    fn the_dosages_of_a_fixed_ploidy_are_the_dosages_of_the_loop_that_reads_the_ploidy() {
        const NUM_INDIVIDUALS: usize = 131;
        for ploidy in 1..=4_usize {
            let of_a_genotype = match NonZeroUsize::new(ploidy) {
                Some(of_a_genotype) => of_a_genotype,
                None => panic!("a ploidy of 0"),
            };
            for (gts, named) in [
                (
                    (0..NUM_INDIVIDUALS * ploidy)
                        .map(|allele| match allele % 11 {
                            0..=3 | 7 | 8 | 10 => 1,
                            4 | 9 => 0,
                            5 => 2,
                            _ => M,
                        })
                        .collect::<Vec<i8>>(),
                    "the genotypes of every dosage",
                ),
                (
                    (0..NUM_INDIVIDUALS * ploidy)
                        .map(|allele| if allele % 11 == 5 { M } else { 1 })
                        .collect::<Vec<i8>>(),
                    "the genotypes of one called dosage",
                ),
            ] {
                let mut of_the_match = (
                    vec![9.0; NUM_INDIVIDUALS],
                    vec![9.0; NUM_INDIVIDUALS],
                    vec![9.0; NUM_INDIVIDUALS],
                );
                let mut of_the_loop = of_the_match.clone();
                let matched = the_dosages_of_a_variant(
                    &gts,
                    of_a_genotype,
                    1,
                    &mut of_the_match.0,
                    &mut of_the_match.1,
                    &mut of_the_match.2,
                );
                let looped = the_dosages_of_any_ploidy(
                    &gts,
                    ploidy,
                    1,
                    &mut of_the_loop.0,
                    &mut of_the_loop.1,
                    &mut of_the_loop.2,
                );
                assert_eq!(
                    of_the_match, of_the_loop,
                    "{named}: the three values of a ploidy of {ploidy}"
                );
                assert_eq!(
                    matched, looped,
                    "{named}: the variance of a ploidy of {ploidy}"
                );
            }
            let gts: Vec<i8> = (0..NUM_INDIVIDUALS * ploidy)
                .map(|allele| match allele % 11 {
                    0..=3 | 7 | 8 | 10 => 1,
                    4 | 9 => 0,
                    5 => 2,
                    _ => M,
                })
                .collect();
            let mut values = (
                vec![9.0; NUM_INDIVIDUALS],
                vec![9.0; NUM_INDIVIDUALS],
                vec![9.0; NUM_INDIVIDUALS],
            );
            let has_variance = the_dosages_of_a_variant(
                &gts,
                of_a_genotype,
                1,
                &mut values.0,
                &mut values.1,
                &mut values.2,
            );
            assert!(
                values.1.contains(&0.0),
                "a genotype with an allele missing at a ploidy of {ploidy}"
            );
            assert!(
                values.0.contains(&0.0) && values.1.contains(&1.0),
                "a genotype of the major allele only at a ploidy of {ploidy}"
            );
            assert!(has_variance, "two dosages at a ploidy of {ploidy}");
            let of_one_dosage: Vec<i8> = (0..NUM_INDIVIDUALS * ploidy)
                .map(|allele| if allele % 11 == 5 { M } else { 1 })
                .collect();
            assert!(
                !the_dosages_of_a_variant(
                    &of_one_dosage,
                    of_a_genotype,
                    1,
                    &mut values.0,
                    &mut values.1,
                    &mut values.2,
                ),
                "one called dosage at a ploidy of {ploidy}"
            );
        }
    }

    #[test]
    fn a_variant_whose_called_genotypes_hold_one_dosage_has_no_variance() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        // v4 has 0/0 in every individual, and the other four have two
        // dosages at least.
        let found: Vec<bool> = (0..5).map(|var| dosages.has_variance(var)).collect();
        assert_eq!(found, vec![true, true, true, false, true]);
        assert!(!dosages.has_variance(5), "a sixth variant has variance");
    }

    #[test]
    fn a_variant_with_no_called_genotype_has_no_dosage_no_variance_and_no_frequency() {
        let block = block_of(&[&[M; 12], &[0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1]], 6, 2);
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(dosages_of(&dosages, 0), vec![None; 6]);
        assert!(!dosages.has_variance(0), "it has variance");
        assert_eq!(dosages.maf(0), None);
        // The variant beside it is v1 of the worked example, so a variant
        // that was called at nothing leaves the next one as it was.
        assert_eq!(
            dosages_of(&dosages, 1),
            vec![Some(0), Some(0), Some(1), Some(1), Some(2), Some(2)]
        );
    }

    #[test]
    fn the_major_allele_frequency_is_the_largest_count_over_the_called_alleles() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        // The counts of the alleles of the five variants of the example:
        // v1 six 0 and six 1, v2 four 0 and eight 1, v3 seven 0 and three
        // 1 of its ten called alleles, v4 twelve 0, and v5 six of each.
        let of_the_counts = [6.0 / 12.0, 8.0 / 12.0, 7.0 / 10.0, 12.0 / 12.0, 6.0 / 12.0];
        for (var, expected) in of_the_counts.iter().enumerate() {
            let found = dosages.maf(var).expect("the frequency");
            assert_the_number_is(
                found,
                *expected,
                &format!("the frequency of the variant {var}"),
            );
        }
        assert_eq!(dosages.maf(5), None, "a sixth variant has a frequency");
    }

    #[test]
    fn the_dosages_of_some_individuals_are_counted_from_the_major_allele_of_those_individuals() {
        let block = the_worked_example();
        // The last two individuals are 1/1 at v1, whose major allele is
        // the 0 over the six, so over these two it is the 1 and both
        // dosages are 0.
        let dosages = LdDosages::of_block(&block, &[4, 5]).expect("the dosages");
        assert_eq!((dosages.num_vars(), dosages.num_individuals()), (5, 2));
        assert_eq!(dosages_of(&dosages, 0), vec![Some(0), Some(0)]);
        assert!(!dosages.has_variance(0), "two genotypes 1/1 have variance");
        assert_the_number_is(
            dosages.maf(0).expect("the frequency"),
            1.0,
            "the frequency of v1 over two individuals of 1/1",
        );
        // The individuals come in the order they were given, and the
        // major allele of v1 over one 1/1 and one 0/0 is the lower
        // numbered of the two that tie.
        let dosages = LdDosages::of_block(&block, &[5, 0]).expect("the dosages");
        assert_eq!(dosages_of(&dosages, 0), vec![Some(2), Some(0)]);
        assert!(dosages.has_variance(0), "0/0 and 1/1 have one dosage");
    }

    #[test]
    fn a_half_called_genotype_counts_its_alleles_and_has_no_dosage() {
        // Four tetraploid individuals: 0/0/0/0, 0/0/1/1, 1/1/1/2 and
        // 0/0/0/., whose called alleles are nine 0, five 1 and one 2, so
        // the major allele is the 0 and the frequency is 9 of 15.
        let block = block_of(&[&[0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 0, 0, 0, M]], 4, 4);
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(
            dosages_of(&dosages, 0),
            vec![Some(0), Some(2), Some(4), None]
        );
        assert!(dosages.has_variance(0), "three dosages are one");
        assert_the_number_is(
            dosages.maf(0).expect("the frequency"),
            9.0 / 15.0,
            "the frequency of a variant with a half called genotype",
        );
    }

    #[test]
    fn the_rows_of_a_range_are_the_dosages_of_those_variants() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let rows = dosages.rows(2, 2).expect("the variants 2 and 3");
        assert_eq!((rows.num_vars(), rows.num_individuals()), (2, 6));
        assert_eq!(
            dosages_of(&rows, 0),
            vec![Some(0), Some(0), Some(0), Some(1), None, Some(2)]
        );
        assert_eq!(dosages_of(&rows, 1), vec![Some(0); 6]);
        assert_eq!(
            (rows.has_variance(0), rows.has_variance(1)),
            (true, false),
            "the variance of v3 and v4"
        );
        assert_the_number_is(
            rows.maf(0).expect("the frequency"),
            7.0 / 10.0,
            "the frequency of v3 among the rows",
        );
        assert_the_values_are(
            rows.called.as_slice(),
            &[1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "the called genotypes of the rows",
        );
        // The last variant, and a range of no variant at the end of them.
        assert_eq!(dosages.rows(4, 1).expect("the last variant").num_vars(), 1);
        assert_eq!(dosages.rows(5, 0).expect("no variant").num_vars(), 0);
    }

    #[test]
    fn the_variants_of_a_range_that_is_not_in_the_dosages_are_refused() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        for (first, asked_for) in [(4, 2), (5, 1), (6, 0), (1, usize::MAX)] {
            match dosages.rows(first, asked_for) {
                Err(Error::LdRowsNotInTheDosages {
                    first: found_first,
                    asked_for: found_asked_for,
                    num_vars,
                }) => {
                    assert_eq!(
                        (found_first, found_asked_for, num_vars),
                        (first, asked_for, 5)
                    );
                }
                other => panic!("the variants {first}..{asked_for} were given: {other:?}"),
            }
        }
    }

    #[test]
    fn a_block_with_variants_and_no_genotypes_is_refused() {
        let mut block = the_worked_example();
        block.gts = Vec::new();
        match LdDosages::of_block(&block, &[]) {
            Err(Error::FieldsNotInTheBlock { fields }) => {
                assert_eq!(fields.to_string(), "`gts`");
            }
            other => panic!("a block of no genotype gave dosages: {other:?}"),
        }
    }

    #[test]
    fn a_block_whose_genotypes_are_not_of_its_size_is_refused() {
        let mut block = the_worked_example();
        block.gts.pop();
        match LdDosages::of_block(&block, &[]) {
            Err(Error::BlockArrayOfAnotherSize {
                array,
                found,
                expected,
            }) => {
                assert_eq!((array, found, expected), ("gts", 59, 60));
            }
            other => panic!("a block of 59 alleles gave dosages: {other:?}"),
        }
    }

    #[test]
    fn an_index_that_is_not_an_individual_of_the_block_is_refused() {
        let block = the_worked_example();
        for asked_for in [6, usize::MAX] {
            match LdDosages::of_block(&block, &[0, asked_for]) {
                Err(Error::LdIndividualNotInTheDataset {
                    individual,
                    num_individuals,
                }) => {
                    assert_eq!((individual, num_individuals), (asked_for, 6));
                }
                other => panic!("the individual {asked_for} of six gave dosages: {other:?}"),
            }
        }
    }

    #[test]
    fn a_genotype_of_more_alleles_than_a_dosage_holds_is_refused() {
        let ploidy = MAX_PLOIDY_OF_THE_DOSAGES
            .checked_add(1)
            .expect("one allele more than a dosage holds");
        let gts = vec![0; ploidy];
        let block = block_of(&[gts.as_slice()], 1, ploidy);
        match LdDosages::of_block(&block, &[]) {
            Err(Error::LdPloidyTooLarge { ploidy: found }) => {
                assert_eq!(found, ploidy);
                let message = Error::LdPloidyTooLarge { ploidy: found }.to_string();
                assert!(message.contains("256 alleles"), "{message}");
                assert!(message.contains("255"), "{message}");
            }
            other => panic!("a ploidy of {ploidy} gave dosages: {other:?}"),
        }
        // The largest ploidy a dosage holds is read, and its genotype of
        // 255 alleles that are not the major one has the dosage 255.
        let mut gts = vec![0; MAX_PLOIDY_OF_THE_DOSAGES];
        gts.extend(vec![1; MAX_PLOIDY_OF_THE_DOSAGES]);
        let block = block_of(&[gts.as_slice()], 2, MAX_PLOIDY_OF_THE_DOSAGES);
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(dosages_of(&dosages, 0), vec![Some(0), Some(255)]);
    }

    #[test]
    fn dosages_of_more_values_than_the_linear_algebra_counts_in_are_refused() {
        // 46340 x 46340 is 2147395600 values, the largest square matrix
        // the linear algebra takes, and 46341 x 46341 is 2147488281,
        // which is above it.
        assert_eq!(the_values_of(46340, 46340), Some(2_147_395_600));
        assert_eq!(the_values_of(46341, 46341), None);
        assert_eq!(
            the_values_of(MAX_VALUES_OF_THE_DOSAGES, 1),
            Some(MAX_VALUES_OF_THE_DOSAGES)
        );
        assert_eq!(the_values_of(MAX_VALUES_OF_THE_DOSAGES, 2), None);
        assert_eq!(the_values_of(usize::MAX, 2), None);
        let message = Error::LdDosagesTooLarge {
            num_vars: 46341,
            num_individuals: 46341,
        }
        .to_string();
        assert!(message.contains("46341 variants"), "{message}");
        assert!(message.contains("2147483647"), "{message}");
    }

    #[test]
    fn a_genotype_that_is_not_in_the_row_is_read_as_missing_and_not_as_the_variant_before_it() {
        // The genotypes of the individuals that were asked for are read
        // into one buffer that is kept from one variant to the next, so a
        // genotype that was not found would otherwise leave the one of the
        // variant before it there and be counted as a called genotype.
        // `of_block` refuses an individual that is not one of the block, so
        // nothing of popnei reaches this.
        let of_a_genotype = NonZeroUsize::new(2).expect("the ploidy");
        let mut chosen = vec![1_i8; 4];
        the_genotypes_of(&[0, 0, 1, 1], &[1, 7], of_a_genotype, &mut chosen);
        assert_eq!(chosen, vec![1, 1, M, M]);
    }

    #[test]
    fn a_block_of_more_alleles_in_one_variant_than_the_r2_comes_out_of_exactly_is_refused() {
        // 94906265 alleles in one variant is sqrt(2^53) rounded down, the
        // largest whose square an f64 holds one by one, and the products
        // the r² takes of its six sums reach that square. At the ploidy of
        // 255 the dosages take it is 372181 individuals, 94906155 alleles,
        // and one individual more is 94906410, which is above it. The
        // blocks are of no variant, so nothing of them is read.
        let refused = block_of(&[], 372_182, MAX_PLOIDY_OF_THE_DOSAGES);
        match LdDosages::of_block(&refused, &[]) {
            Err(Error::LdTooManyAllelesInAVariant {
                num_individuals,
                ploidy,
            }) => {
                assert_eq!((num_individuals, ploidy), (372_182, 255));
                let message = Error::LdTooManyAllelesInAVariant {
                    num_individuals,
                    ploidy,
                }
                .to_string();
                assert!(message.contains("94906265 alleles"), "{message}");
                assert!(message.contains("47453132 diploid"), "{message}");
            }
            other => panic!("94906410 alleles in one variant gave dosages: {other:?}"),
        }
        let taken = block_of(&[], 372_181, MAX_PLOIDY_OF_THE_DOSAGES);
        assert_eq!(
            LdDosages::of_block(&taken, &[])
                .expect("94906155 alleles in one variant")
                .num_individuals(),
            372_181
        );
        // The same bound over the individuals that were asked for: two of
        // a block of that many are taken.
        assert_eq!(
            LdDosages::of_block(&refused, &[0, 1])
                .expect("two individuals of the block")
                .num_individuals(),
            2
        );
        assert_eq!(MAX_ALLELES_OF_A_VARIANT, 94_906_265);
    }

    /// A product of the r² that the linear algebra did not do says which
    /// of the six sums it was and what the linear algebra said.
    ///
    /// Nothing of popnei reaches it: the values of the three matrices are
    /// whole numbers and every dimension is checked where the dosages are
    /// built, so what is left is a result of more values than the routines
    /// of BLAS and LAPACK count in, which is 17 GB of r² that no machine
    /// gives. So what is checked here is the text a user would read.
    #[test]
    fn a_product_the_linear_algebra_refused_says_which_of_the_six_sums_it_was() {
        let message = Error::LdLinalg {
            operation: "Σxy",
            source: LinalgError::Dimension {
                argument: "c",
                expected: "2147483647 values at most".to_owned(),
            },
        }
        .to_string();
        assert!(message.contains("Σxy"), "{message}");
        assert!(message.contains("2147483647 values at most"), "{message}");
        assert!(message.contains("the argument c"), "{message}");
    }

    #[test]
    fn a_matrix_this_machine_has_not_the_memory_for_is_an_error_and_not_the_end_of_the_process() {
        // The memory of every matrix of the r² is asked for with
        // `try_reserve_exact`, which gives it back as an error where
        // `vec![0.0; n]` would end the process, and which refuses a number
        // of values whose bytes this machine does not count before it asks
        // the allocator for anything.
        let values = usize::MAX;
        match a_vector_of(
            0.0_f64,
            values,
            &the_memory_for("the dosages", values, size_of::<f64>()),
        ) {
            Err(Error::LdNoMemory {
                what,
                values: found,
                bytes_per_value,
            }) => {
                assert_eq!(
                    (what, found, bytes_per_value),
                    ("the dosages", values, 8_usize)
                );
            }
            Ok(given) => panic!("{} values of 8 bytes were given", given.len()),
            Err(other) => panic!("the memory failed with another error: {other:?}"),
        }
        let message = Error::LdNoMemory {
            what: "Σxy of the r²",
            values: 25,
            bytes_per_value: 8,
        }
        .to_string();
        assert!(message.contains("Σxy of the r²"), "{message}");
        assert!(message.contains("25 values of 8 bytes"), "{message}");
        // What the pass of the fall-off of r² with distance keeps is not
        // all of it values of 8 bytes: a genotype of the window is one
        // byte, and the message says the size the values have.
        let of_the_genotypes = Error::LdNoMemory {
            what: "the genotypes of the variants of the window",
            values: 1_500_000,
            bytes_per_value: size_of::<i8>(),
        }
        .to_string();
        assert!(
            of_the_genotypes.contains("1500000 values of 1 bytes"),
            "{of_the_genotypes}"
        );
    }

    #[test]
    fn a_block_of_no_variant_gives_dosages_of_no_variant() {
        let block = block_of(&[], 6, 2);
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!((dosages.num_vars(), dosages.num_individuals()), (0, 6));
        assert!(dosages.dosages(0).is_none(), "a variant has dosages");
        assert_eq!(dosages.maf(0), None);
    }

    /// Six of the seven pairs of the worked example of "How it is
    /// verified" of `docs/specs/ld.md`: the two variants of the pair,
    /// counted from 0, its six sums, n, Σx, Σy, Σxy, Σxx and Σyy, and its
    /// r². The seventh, v1 against v4, has no sums in the spec, since v4
    /// has one dosage in every individual, and it is asserted as NaN with
    /// the other pairs of v4.
    const THE_PAIRS_OF_THE_EXAMPLE: [(usize, usize, [f64; 6], f64); 6] = [
        (0, 1, [6.0, 6.0, 4.0, 1.0, 10.0, 6.0], 0.675),
        (0, 2, [5.0, 4.0, 3.0, 5.0, 6.0, 5.0], 0.7544642857142857),
        (0, 4, [6.0, 6.0, 6.0, 5.0, 10.0, 10.0], 0.0625),
        (1, 2, [5.0, 4.0, 3.0, 0.0, 6.0, 5.0], 0.6428571428571429),
        (1, 4, [6.0, 4.0, 6.0, 4.0, 6.0, 10.0], 0.0),
        (2, 4, [5.0, 3.0, 4.0, 1.0, 5.0, 6.0], 0.21875),
    ];

    /// The r² of every pair of the two sets of dosages, the variants of
    /// the first one row after another.
    fn the_r2_of(a: &LdDosages, b: &LdDosages) -> Vec<f64> {
        let num_values = a
            .num_vars()
            .checked_mul(b.num_vars())
            .expect("one value for each pair");
        let mut r2 = vec![0.0; num_values];
        r2_between(a, b, &mut r2).expect("the r²");
        r2
    }

    /// The value of the pair of the variant `of_a` of the first set and
    /// the variant `of_b` of the second, in a matrix that holds
    /// `num_vars_of_b` values in each row.
    fn of_the_pair(matrix: &[f64], num_vars_of_b: usize, of_a: usize, of_b: usize) -> f64 {
        matrix
            .chunks_exact(num_vars_of_b)
            .nth(of_a)
            .and_then(|row| row.get(of_b))
            .copied()
            .unwrap_or_else(|| panic!("the pair {of_a}, {of_b} is not in the matrix"))
    }

    /// That the r² is the one of the spec, within the 1e-12 relative of
    /// "How it is verified", which for the 0 of a pair whose dosages do
    /// not vary together at all is equality.
    fn assert_the_r2_is(found: f64, expected: f64, what: &str) {
        assert!(
            (found - expected).abs() <= 1e-12 * expected.abs(),
            "{what}: the r² is {found} and not {expected}"
        );
    }

    /// The six sums of the pair of the variant `of_a` of the first set and
    /// the variant `of_b` of the second, n, Σx, Σy, Σxy, Σxx and Σyy, in
    /// the order of the table of the spec.
    ///
    /// The sums have to be those of two sets that are not one set against
    /// itself, which `the_six_sums_of` builds: a set against itself takes
    /// Σy and Σyy from Σx and Σxx read the other way round and holds
    /// neither of its own, and what it gives is the r² the tests of the
    /// worked example read.
    fn the_sums_of_the_pair(
        sums: &TheSumsOfThePairs,
        num_vars_of_b: usize,
        of_a: usize,
        of_b: usize,
    ) -> [f64; 6] {
        let TheSumsOfTheSecondSet::OfTheirOwn {
            of_b: of_the_second,
            squares_of_b: squares_of_the_second,
        } = &sums.of_the_second_set
        else {
            panic!("the sums are of one set against itself and hold no Σy of their own")
        };
        [
            of_the_pair(&sums.num_individuals, num_vars_of_b, of_a, of_b),
            of_the_pair(&sums.of_a, num_vars_of_b, of_a, of_b),
            of_the_pair(of_the_second, num_vars_of_b, of_a, of_b),
            of_the_pair(&sums.products, num_vars_of_b, of_a, of_b),
            of_the_pair(&sums.squares_of_a, num_vars_of_b, of_a, of_b),
            of_the_pair(squares_of_the_second, num_vars_of_b, of_a, of_b),
        ]
    }

    /// The six sums of every pair of the variants of `dosages` with
    /// themselves, taken over two sets that hold the same variants and are
    /// not one reference, so that Σy and Σyy are two products of their own
    /// and can be read pair by pair.
    fn the_six_sums_of(dosages: &LdDosages) -> TheSumsOfThePairs {
        let num_vars = dosages.num_vars();
        let of_the_same_variants = dosages.rows(0, num_vars).expect("the same variants");
        let num_values = num_vars.checked_mul(num_vars).expect("one for each pair");
        TheSumsOfThePairs::of(&of_the_same_variants, dosages, num_values).expect("the sums")
    }

    /// That the six sums of a pair are the whole numbers of the spec.
    #[expect(
        clippy::float_cmp,
        reason = "the six sums are whole numbers below 2^53, which an f64 holds exactly, so one that is not the number of the spec is not a rounding"
    )]
    fn assert_the_sums_are(found: &[f64; 6], expected: &[f64; 6], what: &str) {
        let named = ["n", "Σx", "Σy", "Σxy", "Σxx", "Σyy"];
        for (sum, (found, expected)) in named.iter().zip(found.iter().zip(expected)) {
            assert!(
                *found == *expected,
                "{what}: {sum} is {found} and not {expected}"
            );
        }
    }

    /// That the two cells of a pair hold the same r².
    ///
    /// r² is the same whichever variant of the pair comes first, and the
    /// two cells come out of the same six sums the other way round, so
    /// they are the same value and not two values within a tolerance.
    #[expect(
        clippy::float_cmp,
        reason = "the two cells of a pair are the same products of the same six whole numbers, taken in the other order"
    )]
    fn assert_the_r2_is_the_same(found: f64, back: f64, what: &str) {
        assert!(
            found == back || (found.is_nan() && back.is_nan()),
            "{what}: it is {found} one way round and {back} the other"
        );
    }

    /// That a matrix of r² holds the values expected, `None` for a pair
    /// that has none.
    fn assert_the_matrix_is(found: &[f64], expected: &[Option<f64>], what: &str) {
        assert_eq!(
            found.len(),
            expected.len(),
            "{what}: the values are not as many"
        );
        for (at, (found, expected)) in found.iter().zip(expected).enumerate() {
            match expected {
                None => assert!(
                    found.is_nan(),
                    "{what}: the pair {at} is {found} and not NaN"
                ),
                Some(expected) => {
                    assert_the_r2_is(*found, *expected, &format!("{what}: the pair {at}"));
                }
            }
        }
    }

    #[test]
    fn the_r2_of_the_worked_example_is_the_one_of_the_spec() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let matrix = the_r2_of(&dosages, &dosages);
        for (of_a, of_b, _, expected) in THE_PAIRS_OF_THE_EXAMPLE {
            assert_the_r2_is(
                of_the_pair(&matrix, 5, of_a, of_b),
                expected,
                &format!("the pair of the variants {of_a} and {of_b}"),
            );
            // r² is the same whichever variant of the pair comes first.
            assert_the_r2_is(
                of_the_pair(&matrix, 5, of_b, of_a),
                expected,
                &format!("the pair of the variants {of_b} and {of_a}"),
            );
        }
        // v4, the variant 3, has one dosage in every individual, so it has
        // no r² against any variant, itself among them, and the seventh
        // pair of the spec, v1 against v4, is one of these.
        for var in 0..5 {
            let of_v4 = of_the_pair(&matrix, 5, 3, var);
            assert!(
                of_v4.is_nan(),
                "the pair of v4 and the variant {var}: {of_v4}"
            );
            let against_v4 = of_the_pair(&matrix, 5, var, 3);
            assert!(
                against_v4.is_nan(),
                "the pair of the variant {var} and v4: {against_v4}"
            );
        }
        // A variant that has two dosages among its called genotypes has an
        // r² of 1 against itself.
        for var in [0, 1, 2, 4] {
            assert_the_r2_is(
                of_the_pair(&matrix, 5, var, var),
                1.0,
                &format!("the variant {var} against itself"),
            );
        }
    }

    #[test]
    fn the_six_sums_of_the_worked_example_are_the_whole_numbers_of_the_spec() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let sums = the_six_sums_of(&dosages);
        for (of_a, of_b, expected, _) in THE_PAIRS_OF_THE_EXAMPLE {
            let found = the_sums_of_the_pair(&sums, 5, of_a, of_b);
            assert_the_sums_are(
                &found,
                &expected,
                &format!("the pair of the variants {of_a} and {of_b}"),
            );
            // The same pair the other way round has Σx and Σy, and Σxx and
            // Σyy, the other way round too, which is what lets a set
            // against itself read Σy and Σyy from Σx and Σxx.
            let [n, of_x, of_y, products, squares_of_x, squares_of_y] = expected;
            let back = the_sums_of_the_pair(&sums, 5, of_b, of_a);
            assert_the_sums_are(
                &back,
                &[n, of_y, of_x, products, squares_of_y, squares_of_x],
                &format!("the pair of the variants {of_b} and {of_a}"),
            );
        }
    }

    /// Two sets that hold as many variants and are not the same dosages.
    ///
    /// The six products are taken for them, as for any two sets that are
    /// not one set against itself, and the four of a set against itself
    /// would give the pairs of v1 and v2 against v3 and v4 as
    /// 0.413265306122449, 1.2, 2.938775510204082 and 0.64, two of them
    /// above the 1 that an r² reaches at most: Σy and Σyy are the
    /// transposes of Σx and Σxx only when both sets hold the same
    /// variants, and here they hold different variants of the same
    /// number.
    #[test]
    fn the_r2_of_two_sets_of_as_many_variants_that_are_not_the_same_set_is_that_of_their_pairs() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let first_two = dosages.rows(0, 2).expect("v1 and v2");
        let next_two = dosages.rows(2, 2).expect("v3 and v4");
        assert_the_matrix_is(
            &the_r2_of(&first_two, &next_two),
            &[
                Some(0.7544642857142857),
                None,
                Some(0.6428571428571429),
                None,
            ],
            "v1 and v2 against v3 and v4",
        );
    }

    #[test]
    fn the_r2_of_two_sets_of_different_sizes_is_that_of_each_of_their_pairs() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let first_two = dosages.rows(0, 2).expect("v1 and v2");
        let last_three = dosages.rows(2, 3).expect("v3, v4 and v5");
        // Two sets that are not the same dosages, so the six products are
        // taken: the two rows of v1 and of v2 against the three columns of
        // v3, of v4, which has no r², and of v5.
        assert_the_matrix_is(
            &the_r2_of(&first_two, &last_three),
            &[
                Some(0.7544642857142857),
                None,
                Some(0.0625),
                Some(0.6428571428571429),
                None,
                Some(0.0),
            ],
            "the two variants against the three",
        );
        assert_the_matrix_is(
            &the_r2_of(&last_three, &first_two),
            &[
                Some(0.7544642857142857),
                Some(0.6428571428571429),
                None,
                None,
                Some(0.0625),
                Some(0.0),
            ],
            "the three variants against the two",
        );
        // One variant against the five, which is the shape of pyNei's
        // `test_the_r_matrix_does_not_build_the_square_of_both_sets_together`.
        let v3 = dosages.rows(2, 1).expect("v3");
        assert_the_matrix_is(
            &the_r2_of(&v3, &dosages),
            &[
                Some(0.7544642857142857),
                Some(0.6428571428571429),
                Some(1.0),
                None,
                Some(0.21875),
            ],
            "v3 against the five variants",
        );
    }

    #[test]
    fn a_pair_whose_variants_were_called_in_no_individual_together_has_no_r2() {
        // Two variants of four individuals, each called in the two the
        // other was not: they have two dosages each, and no individual to
        // be counted in a pair.
        let block = block_of(
            &[&[0, 0, 0, 1, M, M, M, M], &[M, M, M, M, 0, 0, 0, 1]],
            4,
            2,
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert!(
            dosages.has_variance(0) && dosages.has_variance(1),
            "a variant of 0/0 and 0/1 has one dosage"
        );
        assert_the_sums_are(
            &the_sums_of_the_pair(&the_six_sums_of(&dosages), 2, 0, 1),
            &[0.0; 6],
            "the pair of two variants with no individual in common",
        );
        assert_the_matrix_is(
            &the_r2_of(&dosages, &dosages),
            &[Some(1.0), None, None, Some(1.0)],
            "two variants called in no individual together",
        );
    }

    #[test]
    fn two_variants_with_variance_have_no_r2_when_the_individuals_of_the_pair_hold_one_dosage() {
        // The first variant is 0/0 0/0 0/1 0/1 and the second 0/0 0/1 ./.
        // ./., so both have two dosages among their called genotypes and
        // the two individuals called at both hold 0 at the first.
        let block = block_of(
            &[&[0, 0, 0, 0, 0, 1, 0, 1], &[0, 0, 0, 1, M, M, M, M]],
            4,
            2,
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert!(
            dosages.has_variance(0) && dosages.has_variance(1),
            "one of the two variants has one dosage"
        );
        assert_the_sums_are(
            &the_sums_of_the_pair(&the_six_sums_of(&dosages), 2, 0, 1),
            &[2.0, 0.0, 1.0, 0.0, 0.0, 1.0],
            "the pair of two individuals of one dosage at the first variant",
        );
        assert_the_matrix_is(
            &the_r2_of(&dosages, &dosages),
            &[Some(1.0), None, None, Some(1.0)],
            "two variants whose shared individuals hold one dosage",
        );
    }

    #[test]
    fn the_r2_of_a_set_of_no_variant_writes_nothing() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let none = dosages.rows(5, 0).expect("no variant");
        assert!(the_r2_of(&none, &dosages).is_empty());
        assert!(the_r2_of(&dosages, &none).is_empty());
        assert!(the_r2_of(&none, &none).is_empty());
    }

    #[test]
    fn dosages_of_no_individual_give_no_r2_for_any_pair() {
        // No block of popnei gives these: one of variants and no
        // individual holds no genotype, which `of_block` refuses. Every
        // pair of them has an n of 0.
        let dosages = LdDosages {
            num_vars: 2,
            individuals: Vec::new(),
            dosages: Vec::new(),
            called: Vec::new(),
            squares: Vec::new(),
            has_variance: vec![false; 2],
            maf: vec![None; 2],
        };
        assert_the_matrix_is(
            &the_r2_of(&dosages, &dosages),
            &[None; 4],
            "two variants of no individual",
        );
    }

    #[test]
    fn an_out_that_does_not_hold_one_value_for_each_pair_is_refused() {
        let dosages = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        let two = dosages.rows(0, 2).expect("v1 and v2");
        for (num_values, num_vars_of_b) in [(24, 5), (26, 5), (0, 5), (9, 2)] {
            let against = match num_vars_of_b {
                2 => &two,
                _ => &dosages,
            };
            let mut r2 = vec![0.0; num_values];
            match r2_between(&dosages, against, &mut r2) {
                Err(Error::LdR2OfAnotherSize {
                    num_values: found,
                    num_vars_of_a,
                    num_vars_of_b: found_of_b,
                }) => {
                    assert_eq!(
                        (found, num_vars_of_a, found_of_b),
                        (num_values, 5, num_vars_of_b)
                    );
                }
                other => panic!("a buffer of {num_values} values was taken: {other:?}"),
            }
        }
    }

    #[test]
    fn dosages_of_another_number_of_individuals_are_refused() {
        let block = the_worked_example();
        let of_six = LdDosages::of_block(&block, &[]).expect("the dosages of the six");
        let of_two = LdDosages::of_block(&block, &[0, 1]).expect("the dosages of two");
        let mut r2 = vec![0.0; 25];
        match r2_between(&of_six, &of_two, &mut r2) {
            Err(Error::LdDosagesOfOtherIndividuals {
                problem: TheIndividualsThatDiffer::NotAsMany { of_a, of_b },
            }) => {
                assert_eq!((of_a, of_b), (6, 2));
            }
            other => panic!("dosages of six individuals against two were taken: {other:?}"),
        }
    }

    #[test]
    fn dosages_of_as_many_individuals_that_are_not_the_same_ones_are_refused() {
        let block = the_worked_example();
        let of_the_first_three = LdDosages::of_block(&block, &[0, 1, 2]).expect("the first three");
        let of_the_last_three = LdDosages::of_block(&block, &[3, 4, 5]).expect("the last three");
        let mut r2 = vec![0.0; 25];
        match r2_between(&of_the_first_three, &of_the_last_three, &mut r2) {
            Err(Error::LdDosagesOfOtherIndividuals {
                problem: TheIndividualsThatDiffer::NotTheSame { at, of_a, of_b },
            }) => {
                assert_eq!((at, of_a, of_b), (0, 0, 3));
            }
            other => panic!("dosages of three individuals against three others: {other:?}"),
        }
        // The two sets differ at their second individual and not at their
        // first, which is the one the error names.
        let of_two = LdDosages::of_block(&block, &[0, 1]).expect("the first two");
        let of_two_others = LdDosages::of_block(&block, &[0, 2]).expect("two others");
        let mut r2 = vec![0.0; 25];
        match r2_between(&of_two, &of_two_others, &mut r2) {
            Err(Error::LdDosagesOfOtherIndividuals {
                problem: TheIndividualsThatDiffer::NotTheSame { at, of_a, of_b },
            }) => {
                assert_eq!((at, of_a, of_b), (1, 1, 2));
            }
            other => panic!("dosages of two individuals against two others: {other:?}"),
        }
    }

    #[test]
    fn the_dosages_of_every_individual_are_those_of_all_of_them_named_one_by_one() {
        let block = the_worked_example();
        let of_the_block = LdDosages::of_block(&block, &[]).expect("the dosages of the six");
        let named = LdDosages::of_block(&block, &[0, 1, 2, 3, 4, 5]).expect("the six by name");
        // The two hold the same individuals, so the r² of one against the
        // other is taken and is the matrix of the worked example.
        assert_the_r2_is(
            of_the_pair(&the_r2_of(&of_the_block, &named), 5, 0, 1),
            0.675,
            "v1 and v2 of the dosages of the six against the six by name",
        );
        // The variants of a range hold the individuals of the set they
        // come from.
        let rows = of_the_block.rows(0, 2).expect("v1 and v2");
        assert_the_r2_is(
            of_the_pair(&the_r2_of(&rows, &named), 5, 0, 1),
            0.675,
            "v1 and v2 of two variants against the six by name",
        );
    }

    #[test]
    fn an_individual_asked_for_more_than_once_is_refused() {
        let block = the_worked_example();
        // A population is a set of individuals, and one counted twice
        // would be counted twice in n, in the major allele frequency and
        // in every sum of every pair.
        for (individuals, twice) in [
            (vec![0, 0], 0),
            (vec![0, 1, 2, 3, 4, 5, 0], 0),
            (vec![5, 4, 3, 4], 4),
        ] {
            match LdDosages::of_block(&block, &individuals) {
                Err(Error::LdIndividualAskedForTwice { individual }) => {
                    assert_eq!(individual, twice);
                }
                other => panic!("the individuals {individuals:?} gave dosages: {other:?}"),
            }
        }
    }

    /// How many variants each of the two reference datasets holds:
    /// `tests/reference/ld/ld.vcf.gz`, 500 variants of 100 diploid
    /// individuals, and `tests/reference/vcf/many.vcf`, 500 of 50. Each is
    /// read as one block, since the dosages of a set of variants are those
    /// of one block.
    const NUM_VARS_OF_A_REFERENCE: usize = 500;

    /// How many variants `tests/reference/ld/example.vcf` holds, the
    /// worked example of "How it is verified" of `docs/specs/ld.md`: 5
    /// variants of 6 diploid individuals.
    const THE_VARS_OF_THE_EXAMPLE: usize = 5;

    /// How many pairs of two different variants 500 variants have, and how
    /// many of those of `ld.vcf.gz` plink2 gives an r² for and how many it
    /// gives NaN, from "How it is verified" of `docs/specs/ld.md`.
    const THE_PAIRS_OF_THE_LD_DATASET: (usize, usize, usize) = (124_750, 93_096, 31_654);

    /// The path of one of the files of `tests/reference/ld/`, the dataset
    /// of "How it is verified" of `docs/specs/ld.md` with what plink2 and
    /// pyNei give for it.
    ///
    /// The reference files live at the root of the repository, beside the
    /// script that writes them again, and not inside this crate. The path
    /// is built from the directory of the manifest, so it holds whether
    /// the tests are run with `cargo test --workspace` or with `cargo test
    /// -p popnei`.
    fn the_reference_path(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/reference/ld")
            .join(name)
    }

    /// `tests/reference/vcf/many.vcf`, the 500 variants of 50 diploid
    /// individuals of `docs/specs/io_vcf.md`.
    fn many_vcf() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/vcf/many.vcf")
    }

    /// One block with the `num_vars` variants of the VCF at `path`, read
    /// as diploid and with the variants that failed their FILTER among
    /// them, which is what plink2 and pyNei were given.
    fn the_whole_of(path: &Path, needs: Needs, num_vars: usize) -> Block {
        let named = || path.display().to_string();
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(num_vars),
        };
        let mut reader = VcfReader::from_path(path, options)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = named()));
        reader.set_needs(needs);
        let block = reader
            .next_block()
            .unwrap_or_else(|error| panic!("{path}: {error}", path = named()))
            .unwrap_or_else(|| panic!("{path}: it has no variant", path = named()));
        assert_eq!(
            block.num_vars,
            num_vars,
            "{path}: the first block is not the whole file",
            path = named()
        );
        assert!(
            matches!(reader.next_block(), Ok(None)),
            "{path}: it has more variants than the block took",
            path = named()
        );
        block
    }

    /// The square matrix of r² that plink2 wrote for the dataset `name`
    /// and the identifiers of its rows in the order the matrix has them.
    ///
    /// `<name>.unphased.vcor2.bin` holds one float64 for each pair, row
    /// after row, in the byte order of the machine that wrote it, which is
    /// the little endian of every machine popnei is built on, and
    /// `<name>.unphased.vcor2.bin.vars` the identifier of each row, one
    /// per line. `tests/reference/ld/run_plink2.sh` writes both again.
    fn the_matrix_of_plink2(name: &str) -> (Vec<String>, Vec<f64>) {
        let path = the_reference_path(&format!("{name}.unphased.vcor2.bin"));
        let bytes = std::fs::read(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        let values = bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|value| f64::from_le_bytes(*value))
            .collect();
        let path = the_reference_path(&format!("{name}.unphased.vcor2.bin.vars"));
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        let rows = text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(String::from)
            .collect();
        (rows, values)
    }

    /// The dosages pyNei's `to_012` gave for `many.vcf`, which
    /// `tests/reference/ld/make_reference.py` stored: one line for each
    /// variant, the dosage of each individual separated by tabs, and -1
    /// for a genotype with an allele missing, which has no dosage.
    fn the_dosages_of_pynei(name: &str) -> Vec<Vec<Option<u8>>> {
        let path = the_reference_path(name);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()));
        text.lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split('\t')
                    .map(|value| match value.trim() {
                        "-1" => None,
                        dosage => Some(dosage.parse().unwrap_or_else(|error| {
                            panic!("{name}: `{dosage}` is not a dosage: {error}")
                        })),
                    })
                    .collect()
            })
            .collect()
    }

    /// Every pair of two different variants of
    /// `tests/reference/ld/ld.vcf.gz` read with the VCF reader, against
    /// the matrix plink2 v2.0.0-a.7.7 wrote for the same file: of its
    /// 124750 pairs the 93096 that have an r² agree within the 1e-12
    /// relative of "How it is verified" of `docs/specs/ld.md`, and the
    /// other 31654 are NaN on both sides.
    ///
    /// The tolerance is there for a version of plink2 that works the
    /// expression out in another order. On 22 September 2026 every one of
    /// the 93096 came out of the six whole numbers with the bits plink2
    /// has, so a difference of 1e-13 here is something to look at and not
    /// the noise the tolerance allows for.
    #[test]
    fn every_pair_of_the_ld_dataset_is_the_r2_plink2_gives() {
        let block = the_whole_of(
            &the_reference_path("ld.vcf.gz"),
            Needs::GTS | Needs::ID,
            NUM_VARS_OF_A_REFERENCE,
        );
        let ids = block.id.clone().expect("the identifiers of the variants");
        let (rows, of_plink2) = the_matrix_of_plink2("ld");
        // plink2 writes the rows of its matrix in an order of its own,
        // which the .vars file beside it carries, and the r² of a pair is
        // read at the row and the column its two variants have there.
        assert_eq!(
            ids, rows,
            "the variants of the VCF are not the rows of the matrix of plink2"
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(
            (dosages.num_vars(), dosages.num_individuals()),
            (NUM_VARS_OF_A_REFERENCE, 100)
        );
        let matrix = the_r2_of(&dosages, &dosages);
        assert_eq!(
            matrix.len(),
            of_plink2.len(),
            "the two matrices are not as large"
        );
        // The pairs of two different variants, each one once: the two
        // variants, popnei's r² and plink2's. The diagonal, where a
        // variant is against itself, is left out.
        let pairs: Vec<(usize, usize, f64, f64)> = matrix
            .as_chunks::<NUM_VARS_OF_A_REFERENCE>()
            .0
            .iter()
            .zip(of_plink2.as_chunks::<NUM_VARS_OF_A_REFERENCE>().0)
            .enumerate()
            .flat_map(|(of_a, (row, row_of_plink2))| {
                row.iter()
                    .zip(row_of_plink2)
                    .enumerate()
                    .filter(move |(of_b, _)| *of_b > of_a)
                    .map(move |(of_b, (found, expected))| (of_a, of_b, *found, *expected))
            })
            .collect();
        let (num_pairs, with_an_r2, with_none) = THE_PAIRS_OF_THE_LD_DATASET;
        assert_eq!(pairs.len(), num_pairs, "the pairs of two variants");
        // The two counts below add up to every pair, so a pair that one of
        // the two libraries gave an r² for and the other did not falls in
        // neither of them and fails here.
        assert_eq!(
            pairs
                .iter()
                .filter(|(_, _, found, expected)| !found.is_nan() && !expected.is_nan())
                .count(),
            with_an_r2,
            "the pairs that have an r² in both matrices"
        );
        assert_eq!(
            pairs
                .iter()
                .filter(|(_, _, found, expected)| found.is_nan() && expected.is_nan())
                .count(),
            with_none,
            "the pairs that are NaN in both matrices"
        );
        for (of_a, of_b, found, expected) in &pairs {
            if expected.is_nan() {
                continue;
            }
            assert_the_r2_is(
                *found,
                *expected,
                &format!("the pair of the variants {of_a} and {of_b}"),
            );
            // The pairs above are the upper half of the matrix, so what
            // ties the lower half to them is that the two cells of a pair
            // hold the same r².
            assert_the_r2_is_the_same(
                *found,
                of_the_pair(&matrix, NUM_VARS_OF_A_REFERENCE, *of_b, *of_a),
                &format!("the pair of the variants {of_a} and {of_b}"),
            );
        }
        // The 68 variants of the dataset with no variance have NaN in
        // their row, their column and their diagonal cell, which "Missing
        // genotypes" of `docs/specs/ld.md` asks of this test, and each of
        // the other 432 has an r² of exactly 1 against itself.
        let rows = matrix.as_chunks::<NUM_VARS_OF_A_REFERENCE>().0;
        let mut with_variance = 0_usize;
        let mut without_variance = 0_usize;
        for (var, row) in rows.iter().enumerate() {
            let diagonal = row.get(var).copied().expect("the diagonal cell");
            if dosages.has_variance(var) {
                with_variance = with_variance.checked_add(1).expect("the variants counted");
                assert_the_r2_is(diagonal, 1.0, &format!("the variant {var} against itself"));
                continue;
            }
            without_variance = without_variance
                .checked_add(1)
                .expect("the variants counted");
            assert!(
                row.iter().all(|r2| r2.is_nan()),
                "the variant {var} has no variance and its row is not all NaN"
            );
            assert!(
                rows.iter()
                    .all(|row| row.get(var).is_some_and(|r2| r2.is_nan())),
                "the variant {var} has no variance and its column is not all NaN"
            );
        }
        assert_eq!(
            (with_variance, without_variance),
            (432, 68),
            "the variants of the dataset with variance and without it"
        );
    }

    /// The worked example read from `tests/reference/ld/example.vcf` with
    /// the VCF reader, against the matrix plink2 v2.0.0-a.7.7 wrote for
    /// that same file and against the literals of "How it is verified" of
    /// `docs/specs/ld.md`.
    ///
    /// The other tests of the example build its block by hand from the
    /// table of the spec. This is what ties that table to the file plink2
    /// was run on: the dosages of the block the reader gives are the ones
    /// of the hand built block, and the r² of all 25 cells is plink2's.
    #[test]
    fn the_r2_of_the_example_vcf_is_the_one_plink2_gives_and_the_one_of_the_spec() {
        let block = the_whole_of(
            &the_reference_path("example.vcf"),
            Needs::GTS | Needs::ID,
            THE_VARS_OF_THE_EXAMPLE,
        );
        let ids = block.id.clone().expect("the identifiers of the variants");
        let (rows, of_plink2) = the_matrix_of_plink2("example");
        assert_eq!(
            ids, rows,
            "the variants of the VCF are not the rows of the matrix of plink2"
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(
            (dosages.num_vars(), dosages.num_individuals()),
            (THE_VARS_OF_THE_EXAMPLE, 6)
        );
        // The block the other tests build by hand holds the genotypes of
        // this file.
        let by_hand = LdDosages::of_block(&the_worked_example(), &[]).expect("the dosages");
        for var in 0..THE_VARS_OF_THE_EXAMPLE {
            assert_eq!(
                dosages_of(&dosages, var),
                dosages_of(&by_hand, var),
                "the dosages of the variant {var} of the file and of the block built by hand"
            );
        }
        let matrix = the_r2_of(&dosages, &dosages);
        assert_eq!(
            matrix.len(),
            of_plink2.len(),
            "the two matrices are not as large"
        );
        for (at, (found, expected)) in matrix.iter().zip(&of_plink2).enumerate() {
            match expected.is_nan() {
                true => assert!(found.is_nan(), "the cell {at} is {found} and not NaN"),
                false => assert_the_r2_is(*found, *expected, &format!("the cell {at}")),
            }
        }
        // And the seven pairs of the table of the spec, the six with
        // numbers and the one that holds the variant of one dosage.
        for (of_a, of_b, _, expected) in THE_PAIRS_OF_THE_EXAMPLE {
            assert_the_r2_is(
                of_the_pair(&matrix, THE_VARS_OF_THE_EXAMPLE, of_a, of_b),
                expected,
                &format!("the pair of the variants {of_a} and {of_b} of the file"),
            );
        }
        assert!(
            of_the_pair(&matrix, THE_VARS_OF_THE_EXAMPLE, 0, 3).is_nan(),
            "the pair of v1 and v4 of the file"
        );
    }

    /// The four products of one set of variants against itself and the
    /// six of two sets give the same r², to the bit, over the 250000
    /// pairs of `tests/reference/ld/ld.vcf.gz`.
    ///
    /// The shortcut reads Σy and Σyy of a pair as the Σx and Σxx of the
    /// pair the other way round instead of taking two more products, and
    /// this is what would show the two paths drifting apart. They are the
    /// same number because the sums are whole numbers below 2^53, which
    /// `TheSumsOfTheSecondSet::TheOtherWayRound` says why.
    #[test]
    fn the_four_products_of_a_set_against_itself_give_the_r2_the_six_give() {
        let block = the_whole_of(
            &the_reference_path("ld.vcf.gz"),
            Needs::GTS,
            NUM_VARS_OF_A_REFERENCE,
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        // The same variants as a second set of their own, which is not one
        // reference, so Σy and Σyy are two products of their own.
        let of_the_same_variants = dosages
            .rows(0, NUM_VARS_OF_A_REFERENCE)
            .expect("the same variants");
        let of_four = the_r2_of(&dosages, &dosages);
        let of_six = the_r2_of(&of_the_same_variants, &dosages);
        assert_eq!(
            of_four.len(),
            of_six.len(),
            "the two matrices are not as large"
        );
        for (at, (of_four, of_six)) in of_four.iter().zip(&of_six).enumerate() {
            assert_the_r2_is_the_same(*of_four, *of_six, &format!("the pair {at}"));
        }
    }

    /// The dosages popnei reads from `tests/reference/vcf/many.vcf`
    /// against the ones pyNei's `to_012` gave for it, which
    /// `tests/reference/ld/make_reference.py` stored: every one of the
    /// 25000 genotypes of its 500 variants of 50 individuals.
    ///
    /// It is the check of "How it is verified" of `docs/specs/ld.md` for
    /// the variants of more than two alleles and the half called
    /// genotypes, which plink2 cannot make: plink2 counts no allele of a
    /// half called genotype and popnei counts the called one, as
    /// `docs/specs/pca.md` and pyNei do, so the two pick a different major
    /// allele in some variants of more than two alleles and read other
    /// dosages there.
    #[test]
    fn the_dosages_of_many_vcf_are_the_ones_pynei_gives() {
        let block = the_whole_of(&many_vcf(), Needs::GTS, NUM_VARS_OF_A_REFERENCE);
        // The file is the one the spec describes, so the comparison runs
        // over the genotypes that the rule for the major allele is about:
        // 54 of its 500 variants hold more than two alleles, and 257 of
        // its 25000 genotypes have one allele called and one missing.
        let alleles_per_var = block.alleles_per_var().expect("the alleles of a variant");
        let of_more_than_two_alleles = block
            .gts
            .chunks_exact(alleles_per_var)
            .filter(|variant| {
                let mut alleles: Vec<i8> = variant
                    .iter()
                    .copied()
                    .filter(|allele| *allele != MISSING_ALLELE)
                    .collect();
                alleles.sort_unstable();
                alleles.dedup();
                alleles.len() > 2
            })
            .count();
        assert_eq!(
            of_more_than_two_alleles, 54,
            "the variants of more than two alleles"
        );
        let half_called = block
            .gts
            .as_chunks::<2>()
            .0
            .iter()
            .filter(|genotype| {
                genotype
                    .iter()
                    .filter(|allele| **allele == MISSING_ALLELE)
                    .count()
                    == 1
            })
            .count();
        assert_eq!(half_called, 257, "the half called genotypes");
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!(
            (dosages.num_vars(), dosages.num_individuals()),
            (NUM_VARS_OF_A_REFERENCE, 50)
        );
        let of_pynei = the_dosages_of_pynei("many.pynei.dosages.tsv");
        assert_eq!(
            of_pynei.len(),
            NUM_VARS_OF_A_REFERENCE,
            "the variants pyNei read"
        );
        assert_eq!(
            of_pynei.iter().flatten().count(),
            25_000,
            "the genotypes pyNei read"
        );
        for (var, expected) in of_pynei.iter().enumerate() {
            assert_eq!(dosages_of(&dosages, var), *expected, "the variant {var}");
        }
    }

    /// The five pairs of the table of "How it is verified" of
    /// `docs/specs/ld.md`, which plink2 v2.0.0-a.7.7 gave for
    /// `tests/reference/ld/ld.vcf.gz` on 22 September 2026: the
    /// chromosome and the position of each of the two variants, their r²
    /// and the individuals both of them were called in.
    ///
    /// The last pair is on two chromosomes, which `calc_r2_matrix` gives
    /// like any other.
    const THE_PAIRS_OF_THE_TABLE: [(&str, u64, &str, u64, f64, f64); 5] = [
        ("chr1", 1000, "chr1", 2000, 0.353_466_669_239_891, 94.0),
        ("chr1", 1000, "chr1", 3000, 0.398_499_910_809_105_63, 94.0),
        ("chr1", 1000, "chr1", 11000, 0.240_537_842_616_089_57, 93.0),
        ("chr1", 1000, "chr1", 250_000, 0.025_675_192_781_022_8, 95.0),
        ("chr1", 1000, "chr2", 1000, 0.008_140_034_754_693_937, 95.0),
    ];

    /// How many variants of `tests/reference/ld/ld.vcf.gz` have no
    /// variance and how many have some, from "How it is verified" of
    /// `docs/specs/ld.md`.
    const THE_VARIANTS_OF_THE_LD_DATASET: (usize, usize) = (432, 68);

    /// A reader over `tests/reference/ld/ld.vcf.gz`, read as diploid and
    /// with the variants that failed their FILTER among them, which is
    /// what plink2 was given, in blocks of `num_vars_per_block` variants.
    fn the_ld_dataset(
        num_vars_per_block: Option<usize>,
    ) -> VcfReader<std::io::BufReader<std::fs::File>> {
        let path = the_reference_path("ld.vcf.gz");
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block,
        };
        VcfReader::from_path(&path, options)
            .unwrap_or_else(|error| panic!("{path}: {error}", path = path.display()))
    }

    /// The variant of the matrix that is at `pos` of the chromosome
    /// `chrom`.
    fn the_variant_at(matrix: &R2Matrix, chrom: &str, pos: u64) -> usize {
        let names = u32::try_from(matrix.chrom_table().len()).expect("the chromosomes");
        let number = (0..names)
            .find(|number| matrix.chrom_table().name(*number) == Some(chrom))
            .unwrap_or_else(|| panic!("the matrix has no chromosome {chrom}"));
        matrix
            .chroms()
            .iter()
            .zip(matrix.poss())
            .position(|(of_the_var, at)| *of_the_var == number && *at == pos)
            .unwrap_or_else(|| panic!("the matrix has no variant at {chrom}:{pos}"))
    }

    /// The r² of the pair of the variants `of_a` and `of_b` of the matrix.
    fn the_r2_of_the_pair(matrix: &R2Matrix, of_a: usize, of_b: usize) -> f64 {
        of_the_pair(matrix.r2(), matrix.num_vars(), of_a, of_b)
    }

    /// How many individuals both variants of the pair were called in,
    /// which is the n of the table of the spec: the first of the six sums
    /// of the pair, taken over two sets of one variant each that are not
    /// one set against itself.
    fn the_individuals_of_the_pair(dosages: &LdDosages, of_a: usize, of_b: usize) -> f64 {
        let of_a = dosages.rows(of_a, 1).expect("the variant");
        let of_b = dosages.rows(of_b, 1).expect("the variant");
        let sums = TheSumsOfThePairs::of(&of_a, &of_b, 1).expect("the sums");
        sums.num_individuals
            .first()
            .copied()
            .expect("the individuals of the pair")
    }

    /// That two matrices of r² hold the same values, to the bit, with NaN
    /// where both have none.
    #[expect(
        clippy::float_cmp,
        reason = "the two matrices are the same six whole numbers put through the same operations, so a value that is not the same value is not a rounding"
    )]
    fn assert_the_matrices_are_the_same(found: &R2Matrix, expected: &R2Matrix, what: &str) {
        assert_eq!(
            found.num_vars(),
            expected.num_vars(),
            "{what}: the variants are not as many"
        );
        assert_eq!(found.chroms(), expected.chroms(), "{what}: the chromosomes");
        assert_eq!(found.poss(), expected.poss(), "{what}: the positions");
        for (at, (found, expected)) in found.r2().iter().zip(expected.r2()).enumerate() {
            assert!(
                found == expected || (found.is_nan() && expected.is_nan()),
                "{what}: the pair {at} is {found} and not {expected}"
            );
        }
    }

    /// The matrix of every pair of `tests/reference/ld/ld.vcf.gz` read
    /// with the VCF reader: the five pairs of the table of "How it is
    /// verified" of `docs/specs/ld.md` with their n, every one of its
    /// 250000 cells against the matrix plink2 v2.0.0-a.7.7 wrote for the
    /// same file, and the 68 variants of one dosage, whose row, whose
    /// column and whose diagonal cell are NaN.
    ///
    /// The rows of plink2's matrix are the variants of the file in the
    /// order it holds them, which
    /// `every_pair_of_the_ld_dataset_is_the_r2_plink2_gives` asserts
    /// against their identifiers, and `calc_r2_matrix` gives the variants
    /// in the order its reader gave them.
    #[test]
    fn the_matrix_of_the_ld_dataset_is_the_one_plink2_gives() {
        let mut reader = the_ld_dataset(None);
        let matrix =
            calc_r2_matrix(&mut reader, MAX_NUM_VARS_OF_THE_MATRIX).expect("the matrix of r²");
        assert_eq!(matrix.num_vars(), NUM_VARS_OF_A_REFERENCE);
        assert_eq!(matrix.r2().len(), 250_000, "the cells of the matrix");
        assert_eq!(matrix.chroms().len(), NUM_VARS_OF_A_REFERENCE);
        assert_eq!(matrix.poss().len(), NUM_VARS_OF_A_REFERENCE);
        assert_eq!(matrix.chrom_table().len(), 2, "the chromosomes of the file");
        // The five pairs of the table, with the n of each: the r² is read
        // off the matrix and the n off the six sums of the two variants.
        let block = the_whole_of(
            &the_reference_path("ld.vcf.gz"),
            Needs::GTS,
            NUM_VARS_OF_A_REFERENCE,
        );
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        for (chrom_of_a, of_a, chrom_of_b, of_b, r2, individuals) in THE_PAIRS_OF_THE_TABLE {
            let named = format!("the pair of {chrom_of_a}:{of_a} and {chrom_of_b}:{of_b}");
            let of_a = the_variant_at(&matrix, chrom_of_a, of_a);
            let of_b = the_variant_at(&matrix, chrom_of_b, of_b);
            assert_the_r2_is(the_r2_of_the_pair(&matrix, of_a, of_b), r2, &named);
            assert_the_r2_is_the_same(
                the_r2_of_the_pair(&matrix, of_a, of_b),
                the_r2_of_the_pair(&matrix, of_b, of_a),
                &named,
            );
            let found = the_individuals_of_the_pair(&dosages, of_a, of_b);
            assert!(
                (found - individuals).abs() < 0.5,
                "{named}: n is {found} and not {individuals}"
            );
        }
        // And every cell of the matrix against plink2's, the diagonal
        // among them.
        let (_, of_plink2) = the_matrix_of_plink2("ld");
        assert_eq!(
            matrix.r2().len(),
            of_plink2.len(),
            "the two matrices are not as large"
        );
        for (at, (found, expected)) in matrix.r2().iter().zip(&of_plink2).enumerate() {
            match expected.is_nan() {
                true => assert!(found.is_nan(), "the cell {at} is {found} and not NaN"),
                false => assert_the_r2_is(*found, *expected, &format!("the cell {at}")),
            }
        }
        // The variants of one dosage have NaN in their row, their column
        // and their diagonal cell, and each of the others has an r² of 1
        // against itself.
        let rows = matrix.r2().as_chunks::<NUM_VARS_OF_A_REFERENCE>().0;
        let mut with_variance = 0_usize;
        let mut without_variance = 0_usize;
        for (var, row) in rows.iter().enumerate() {
            let diagonal = row.get(var).copied().expect("the diagonal cell");
            if row.iter().all(|r2| r2.is_nan()) {
                without_variance = without_variance.checked_add(1).expect("the variants");
                assert!(
                    rows.iter()
                        .all(|row| row.get(var).is_some_and(|r2| r2.is_nan())),
                    "the row of the variant {var} is all NaN and its column is not"
                );
                continue;
            }
            with_variance = with_variance.checked_add(1).expect("the variants");
            assert_the_r2_is(diagonal, 1.0, &format!("the variant {var} against itself"));
        }
        assert_eq!(
            (with_variance, without_variance),
            THE_VARIANTS_OF_THE_LD_DATASET,
            "the variants of the dataset with variance and without it"
        );
    }

    /// The matrix of `tests/reference/ld/ld.vcf.gz` is the same, to the
    /// bit, for blocks of 7, 64, 256 and 500 variants and for tiles of
    /// those sizes, which is 16 ways of cutting the same 500 variants.
    ///
    /// The blocks are what a user chooses and the tiles are popnei's own.
    /// Neither changes a sum: the six sums of a pair are whole numbers
    /// that an `f64` holds exactly and each of them runs over the
    /// individuals, which no block and no tile cuts.
    #[test]
    fn neither_the_blocks_nor_the_tiles_change_the_matrix() {
        // `THE_VARS_OF_A_TILE` is in the list so that the size the
        // calculation runs at by default is one of the sizes compared, and
        // a performance review that moves it cannot move it out of this
        // test. It is larger than the variants of this dataset, so it is
        // also the case of one tile that holds every variant.
        let sizes = [7, 64, 256, NUM_VARS_OF_A_REFERENCE, THE_VARS_OF_A_TILE];
        let mut reader = the_ld_dataset(Some(NUM_VARS_OF_A_REFERENCE));
        let of_one_block = the_r2_matrix_in_tiles_of(&mut reader, 5000, THE_VARS_OF_A_TILE)
            .expect("the matrix of r²");
        for num_vars_per_block in sizes {
            for vars_per_tile in sizes {
                let mut reader = the_ld_dataset(Some(num_vars_per_block));
                let matrix = the_r2_matrix_in_tiles_of(&mut reader, 5000, vars_per_tile)
                    .expect("the matrix of r²");
                assert_the_matrices_are_the_same(
                    &matrix,
                    &of_one_block,
                    &format!("blocks of {num_vars_per_block} and tiles of {vars_per_tile}"),
                );
            }
        }
    }

    /// The products run on the threads of the backend of the linear
    /// algebra, which is called from outside rayon, and the matrix is the
    /// same to the bit on a pool of one thread and on one of four.
    ///
    /// The pools are built here and are not rayon's global one, which has
    /// one thread per core of the machine. rayon is a dependency of the
    /// targets that are not wasm, so this test is compiled for those
    /// alone.
    #[cfg(not(target_family = "wasm"))]
    #[test]
    fn the_number_of_threads_does_not_change_the_matrix() {
        let in_a_pool = |threads| {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("the pool");
            pool.install(|| {
                let mut reader = the_ld_dataset(Some(64));
                calc_r2_matrix(&mut reader, MAX_NUM_VARS_OF_THE_MATRIX).expect("the matrix of r²")
            })
        };

        let on_one = in_a_pool(1);
        assert_eq!(on_one.num_vars(), NUM_VARS_OF_A_REFERENCE);
        assert_the_matrices_are_the_same(&in_a_pool(4), &on_one, "four threads against one");
    }

    /// The worked example of "How it is verified" of `docs/specs/ld.md`,
    /// given as two blocks of three and two variants and taken in tiles of
    /// two, so that the variants of a tile come from two blocks and the
    /// last tile holds one variant.
    #[test]
    fn the_matrix_of_the_worked_example_is_the_one_of_the_spec() {
        let rows: Vec<&[i8]> = THE_WORKED_EXAMPLE
            .iter()
            .map(|row| row.as_slice())
            .collect();
        let mut first = block_of(&rows[..3], 6, 2);
        let mut second = block_of(&rows[3..], 6, 2);
        the_chrom_and_the_pos_of(&mut first, 0);
        the_chrom_and_the_pos_of(&mut second, 3);
        let mut reader = GivenBlocks::of(vec![first, second], 6, 2);

        let matrix = the_r2_matrix_in_tiles_of(&mut reader, 5000, 2).expect("the matrix of r²");

        assert_eq!(matrix.num_vars(), THE_VARS_OF_THE_EXAMPLE);
        assert_eq!(matrix.chroms(), [0, 0, 0, 0, 0]);
        assert_eq!(matrix.poss(), [1000, 2000, 3000, 4000, 5000]);
        assert_eq!(matrix.chrom_table().name(0), Some("chr1"));
        for (of_a, of_b, _, expected) in THE_PAIRS_OF_THE_EXAMPLE {
            let named = format!("the pair of the variants {of_a} and {of_b}");
            assert_the_r2_is(the_r2_of_the_pair(&matrix, of_a, of_b), expected, &named);
            assert_the_r2_is_the_same(
                the_r2_of_the_pair(&matrix, of_a, of_b),
                the_r2_of_the_pair(&matrix, of_b, of_a),
                &named,
            );
        }
        // v4 has one dosage in every individual, so its row, its column
        // and its diagonal cell are NaN, and the other four variants have
        // an r² of 1 against themselves.
        for var in 0..THE_VARS_OF_THE_EXAMPLE {
            let r2 = the_r2_of_the_pair(&matrix, 3, var);
            assert!(r2.is_nan(), "the pair of v4 and the variant {var} is {r2}");
            let r2 = the_r2_of_the_pair(&matrix, var, 3);
            assert!(r2.is_nan(), "the pair of the variant {var} and v4 is {r2}");
            if var != 3 {
                assert_the_r2_is(
                    the_r2_of_the_pair(&matrix, var, var),
                    1.0,
                    &format!("the variant {var} against itself"),
                );
            }
        }
    }

    /// A pass of more variants than the calculation was allowed is
    /// refused, with both numbers and the memory the matrix would have
    /// needed, at the block that passes the number and not at the end of
    /// the source.
    #[test]
    fn more_variants_than_the_calculation_was_allowed_are_refused() {
        let mut reader = the_ld_dataset(Some(64));

        let error = calc_r2_matrix(&mut reader, 100).expect_err("the variants are more than 100");

        // The second block of 64 variants is where the pass passes the
        // 100, and the matrix of 128 variants is 128 x 128 values of 8
        // bytes.
        let message = error.to_string();
        assert!(
            matches!(
                error,
                Error::LdTooManyVars {
                    num_vars: 128,
                    max_num_vars: 100,
                    bytes: 131_072
                }
            ),
            "the variants above the cap gave: {error:?}"
        );
        assert!(
            message.contains("128") && message.contains("100") && message.contains("131072"),
            "the message holds neither both numbers nor the memory: {message}"
        );
    }

    /// A `max_num_vars` whose matrix holds more values than this machine
    /// counts is refused before the source is read.
    #[test]
    fn a_max_num_vars_whose_matrix_is_not_counted_is_refused() {
        let mut reader = GivenBlocks::of(vec![the_worked_example()], 6, 2);

        let error = calc_r2_matrix(&mut reader, usize::MAX).expect_err("the matrix is not counted");

        assert!(
            matches!(error, Error::LdMaxNumVarsTooLarge { max_num_vars } if max_num_vars == usize::MAX),
            "the cap that is not counted gave: {error:?}"
        );
        assert_eq!(reader.calls, 0, "the source was read");
    }

    /// A reader with no variant is an error and not a matrix of no cell.
    /// It is the one case every calculation over a pass raises, which
    /// carries the counts of the filters of the reader; this reader has
    /// none, so the error says that the source of the pass gave no variant.
    #[test]
    fn a_reader_with_no_variant_is_an_error() {
        let mut reader = GivenBlocks::of(Vec::new(), 6, 2);

        let error = calc_r2_matrix(&mut reader, 5000).expect_err("the reader has no variant");

        assert!(
            matches!(
                &error,
                Error::PassGaveNoVariant {
                    num_vars_of_the_source: 0,
                    filters,
                } if filters.is_empty()
            ),
            "the reader with no variant gave: {error:?}"
        );
    }

    /// A block with variants and no position is the error of a field that
    /// is not in the block: the matrix carries the chromosome and the
    /// position of each of its variants.
    #[test]
    fn a_block_with_no_position_is_the_error_of_a_field_that_is_not_there() {
        let mut reader = GivenBlocks::of(vec![the_worked_example()], 6, 2);

        let error = calc_r2_matrix(&mut reader, 5000).expect_err("the block has no position");

        assert!(
            matches!(error, Error::FieldsNotInTheBlock { fields } if fields == Needs::CHROM_POS),
            "the block with no position gave: {error:?}"
        );
    }

    /// The matrix given away holds what the accessors lend, and the values
    /// are the same ones: a binding crate hands the vector to its language
    /// instead of copying 200 MB out of it.
    #[test]
    fn the_matrix_given_away_holds_what_the_accessors_lend_and_copies_nothing() {
        let mut block = the_worked_example();
        the_chrom_and_the_pos_of(&mut block, 0);
        let mut reader = GivenBlocks::of(vec![block], 6, 2);
        let matrix = calc_r2_matrix(&mut reader, 5000).expect("the matrix of r²");
        let num_vars = matrix.num_vars();
        let chroms = matrix.chroms().to_vec();
        let poss = matrix.poss().to_vec();
        let of_the_first_chrom = matrix.chrom_table().name(0).map(str::to_owned);
        let values = matrix.r2().to_vec();
        // Where the values of the matrix are: the vector given away is the
        // one the matrix held and not a copy of it.
        let where_they_are = matrix.r2().as_ptr();

        let given = matrix.given_away();

        assert_eq!(given.num_vars, num_vars);
        // A NaN is no value's equal, and the matrix of the worked example
        // has five of them, so the values are compared by their bits.
        assert_eq!(given.r2.len(), values.len());
        assert!(
            given
                .r2
                .iter()
                .zip(&values)
                .all(|(given, lent)| given.to_bits() == lent.to_bits()),
            "the values given away are not the ones the accessor lent"
        );
        assert!(std::ptr::eq(given.r2.as_ptr(), where_they_are));
        assert_eq!(given.chroms, chroms);
        assert_eq!(given.poss, poss);
        assert_eq!(
            given.chrom_table.name(0).map(str::to_owned),
            of_the_first_chrom
        );
    }

    /// The calculation asks its reader for the genotypes, the chromosome
    /// and the position, and for nothing else, so a reader over a file
    /// leaves the other columns of a variant unparsed.
    #[test]
    fn the_calculation_asks_its_reader_for_the_genotypes_the_chromosome_and_the_position() {
        let mut block = the_worked_example();
        the_chrom_and_the_pos_of(&mut block, 0);
        let mut reader = GivenBlocks::of(vec![block], 6, 2);

        calc_r2_matrix(&mut reader, 5000).expect("the matrix of r²");

        assert_eq!(reader.needs, Needs::GTS | Needs::CHROM_POS);
    }

    /// A reader that gives a block of other individuals than it says its
    /// source has is a reader with a defect, and the rows of its blocks
    /// cannot be put together into the tiles of one pass.
    #[test]
    fn a_block_of_other_individuals_than_the_reader_says_is_an_error() {
        let mut block = the_worked_example();
        the_chrom_and_the_pos_of(&mut block, 0);
        let mut reader = GivenBlocks::of(vec![block], 3, 2);

        let error = calc_r2_matrix(&mut reader, 5000).expect_err("the block is of others");

        assert!(
            matches!(
                error,
                Error::BlocksDoNotFitTogether {
                    num_individuals: 3,
                    found_num_individuals: 6,
                    ..
                }
            ),
            "the block of other individuals gave: {error:?}"
        );
    }

    /// The error of the reader is given on as it is, and the variants it
    /// gave before it are dropped with the calculation.
    #[test]
    fn the_error_of_the_reader_is_given_on() {
        let mut block = the_worked_example();
        the_chrom_and_the_pos_of(&mut block, 0);
        let mut reader = GivenBlocks::failing_at(vec![block], 6, 2);

        let error = calc_r2_matrix(&mut reader, 5000).expect_err("the reader failed");

        assert!(
            error.to_string().contains(THE_READER_FAILED),
            "the error of the reader gave: {error:?}"
        );
    }

    /// The chromosome `chr1` and the positions 1000, 2000 and on, of the
    /// variant `first` of the dataset up, written into the block.
    fn the_chrom_and_the_pos_of(block: &mut Block, first: u64) {
        let num_vars = u64::try_from(block.num_vars).expect("the variants");
        block.chrom = Some(vec![0; block.num_vars]);
        block.pos = Some(
            (first..first.saturating_add(num_vars))
                .map(|var| var.saturating_add(1).saturating_mul(1000))
                .collect(),
        );
    }

    /// What the reader of these tests says when a test asked it to fail.
    const THE_READER_FAILED: &str = "the reader of the tests failed";

    /// A reader of blocks written for these tests: it gives the blocks it
    /// was built with, says how many individuals and what ploidy its
    /// source has, keeps what it was last asked to fill, and gives an
    /// error instead of its first block when a test asks for one.
    struct GivenBlocks {
        individuals: Vec<String>,
        ploidy: usize,
        chroms: ChromTable,
        /// The blocks still to give, the last one first.
        left: Vec<Block>,
        /// Whether it gives an error instead of a block.
        fails: bool,
        /// How many times it was asked for a block.
        calls: usize,
        /// What it was last asked to fill.
        needs: Needs,
    }

    impl GivenBlocks {
        /// A reader of `num_individuals` individuals of the ploidy
        /// `ploidy`, named `i000` and on, that gives `blocks` in their
        /// order and holds the one chromosome `chr1`.
        fn of(blocks: Vec<Block>, num_individuals: usize, ploidy: usize) -> GivenBlocks {
            let mut left = blocks;
            left.reverse();
            let mut chroms = ChromTable::new();
            chroms.intern("chr1");
            GivenBlocks {
                individuals: (0..num_individuals).map(|at| format!("i{at:03}")).collect(),
                ploidy,
                chroms,
                left,
                fails: false,
                calls: 0,
                needs: Needs::ALL,
            }
        }

        /// The same reader, whose first call is an error.
        fn failing_at(blocks: Vec<Block>, num_individuals: usize, ploidy: usize) -> GivenBlocks {
            GivenBlocks {
                fails: true,
                ..GivenBlocks::of(blocks, num_individuals, ploidy)
            }
        }
    }

    impl BlockReader for GivenBlocks {
        fn next_block(&mut self) -> crate::error::Result<Option<Block>> {
            self.calls = self.calls.saturating_add(1);
            if self.fails {
                return Err(Error::Io(std::io::Error::other(THE_READER_FAILED)));
            }
            Ok(self.left.pop())
        }

        fn individuals(&self) -> &[String] {
            &self.individuals
        }

        fn ploidy(&self) -> usize {
            self.ploidy
        }

        fn chroms(&self) -> &ChromTable {
            &self.chroms
        }

        fn set_needs(&mut self, needs: Needs) {
            self.needs = needs;
        }

        fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
            Vec::new()
        }
    }
}
