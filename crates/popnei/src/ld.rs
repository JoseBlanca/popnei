//! Linkage disequilibrium: how much the genotype of one variant says about
//! the genotype of another.
//!
//! Two variants are in linkage disequilibrium when they sit close enough on
//! a chromosome that few recombinations have separated them, and the
//! measure of it is r², the square of the correlation between the dosages
//! of the two variants over the individuals called at both.
//! `docs/specs/ld.md` is the spec of the module.
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

use popnei_linalg::product;

use crate::block::Block;
use crate::error::{Error, Result};
use crate::variant::{AlleleCounts, MISSING_ALLELE, Needs, count_alleles, the_major_allele};

/// The most values one of the matrices of [`LdDosages`] holds, the variants
/// it was built over times the individuals.
///
/// It is what the routines of BLAS and LAPACK count the values of a matrix
/// in, and `crates/popnei-linalg` refuses a matrix above it on both of its
/// backends, so dosages that no product could be taken over are refused
/// where they are built and not at the first product.
pub const MAX_VALUES_OF_THE_DOSAGES: usize = 2_147_483_647;

/// The most alleles a genotype of [`LdDosages`] holds.
///
/// A dosage is how many alleles of a genotype are not the major allele of
/// its variant, so it is at most the ploidy, and [`LdDosages::dosages`]
/// gives it as one byte. The VCF reader takes 255 alleles in a genotype at
/// most, and the largest ploidy of an organism is a dozen, so this refuses
/// no dataset that a reader of popnei gives.
pub const MAX_PLOIDY_OF_THE_DOSAGES: usize = 255;

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
    /// missing written as 0. It is the A of the spec.
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
    /// individual asked for more than once, a block whose genotypes hold
    /// more than
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
        let mut asked_for_already =
            a_vector_of(false, block.num_individuals, &|| Error::LdNoMemory {
                what: "the individuals asked for",
                values: block.num_individuals,
            })?;
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
        let num_individuals = match individuals.is_empty() {
            true => block.num_individuals,
            false => individuals.len(),
        };
        let too_large = || Error::LdDosagesTooLarge {
            num_vars: block.num_vars,
            num_individuals,
        };
        let values = the_values_of(block.num_vars, num_individuals).ok_or_else(too_large)?;
        let mut dosages = LdDosages {
            num_vars: block.num_vars,
            individuals: match individuals.is_empty() {
                true => (0..block.num_individuals).collect(),
                false => individuals.to_vec(),
            },
            dosages: a_vector_of(0.0, values, &the_memory_for("the dosages", values))?,
            called: a_vector_of(0.0, values, &the_memory_for("the called genotypes", values))?,
            squares: a_vector_of(
                0.0,
                values,
                &the_memory_for("the squares of the dosages", values),
            )?,
            has_variance: vec![false; block.num_vars],
            maf: vec![None; block.num_vars],
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
            false => vec![
                MISSING_ALLELE;
                num_individuals
                    .checked_mul(of_a_genotype.get())
                    .ok_or_else(too_large)?
            ],
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
    /// numbers and how many variants there are.
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
        Ok(LdDosages {
            num_vars,
            individuals: self.individuals.clone(),
            dosages: the_values_from(&self.dosages, from, values),
            called: the_values_from(&self.called, from, values),
            squares: the_values_from(&self.squares, from, values),
            has_variance: self
                .has_variance
                .iter()
                .skip(first)
                .take(num_vars)
                .copied()
                .collect(),
            maf: self
                .maf
                .iter()
                .skip(first)
                .take(num_vars)
                .copied()
                .collect(),
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
/// the matrix of r² against itself is asked for: n and Σxy are then each
/// their own transpose, and Σy and Σyy are the transposes of Σx and Σxx.
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
    let values = out
        .iter_mut()
        .zip(&sums.num_individuals)
        .zip(&sums.products)
        .zip(&sums.of_a)
        .zip(&sums.of_b)
        .zip(&sums.squares_of_a)
        .zip(&sums.squares_of_b);
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
    Ok(())
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
    /// Σy, the sum of the dosages of the variant of the second set over
    /// them.
    of_b: Vec<f64>,
    /// Σxx, the sum of the squares of the dosages of the variant of the
    /// first set over them.
    squares_of_a: Vec<f64>,
    /// Σyy, the sum of the squares of the dosages of the variant of the
    /// second set over them.
    squares_of_b: Vec<f64>,
}

impl TheSumsOfThePairs {
    /// The six sums of every pair of `a` and `b`, which hold the same
    /// individuals in the same order and have `num_values` pairs between
    /// them, the variants of `a` times those of `b`, a number the caller
    /// has counted.
    ///
    /// A product sums over the columns of its first matrix and the rows of
    /// its second, and the three matrices of a set of dosages are variants
    /// x individuals, so the matrices of `b` are transposed to individuals
    /// x variants before they are multiplied. Each transpose copies 8
    /// bytes for every variant of `b` and individual, 4.1 MB for 512
    /// variants of 1000 individuals.
    ///
    /// # Errors
    ///
    /// [`Error::LdLinalg`] when a product could not be worked out, and
    /// [`Error::LdNoMemory`] when this machine did not give the memory of
    /// one of the six sums or of a transpose.
    fn of(a: &LdDosages, b: &LdDosages, num_values: usize) -> Result<TheSumsOfThePairs> {
        let (rows, inner, cols) = (a.num_vars, a.num_individuals(), b.num_vars);
        let called_of_b = the_transpose_of(
            &b.called,
            cols,
            inner,
            "the transpose of the called genotypes of b",
        )?;
        let dosages_of_b =
            the_transpose_of(&b.dosages, cols, inner, "the transpose of the dosages of b")?;
        let mut num_individuals = a_vector_of(0.0, num_values, &the_memory_for("n", num_values))?;
        let mut products = a_vector_of(0.0, num_values, &the_memory_for("Σxy", num_values))?;
        let mut of_a = a_vector_of(0.0, num_values, &the_memory_for("Σx", num_values))?;
        let mut squares_of_a = a_vector_of(0.0, num_values, &the_memory_for("Σxx", num_values))?;
        let sum_of = |of_the_variants: &[f64], by_individual: &[f64], into: &mut [f64], sum| {
            product(of_the_variants, rows, inner, by_individual, cols, into).map_err(|source| {
                Error::LdLinalg {
                    operation: sum,
                    source,
                }
            })
        };
        sum_of(&a.called, &called_of_b, &mut num_individuals, "n")?;
        sum_of(&a.dosages, &dosages_of_b, &mut products, "Σxy")?;
        sum_of(&a.dosages, &called_of_b, &mut of_a, "Σx")?;
        sum_of(&a.squares, &called_of_b, &mut squares_of_a, "Σxx")?;
        let (of_b, squares_of_b) = if std::ptr::eq(a, b) {
            // One set of variants against itself: the pair of the variants
            // i and j holds the two variants of the pair of j and i the
            // other way round, so Σy and Σyy are the transposes of Σx and
            // Σxx and two of the six products are not taken.
            (
                the_transpose_of(&of_a, rows, cols, "the transpose of Σx")?,
                the_transpose_of(&squares_of_a, rows, cols, "the transpose of Σxx")?,
            )
        } else {
            let squares_of_b_by_individual = the_transpose_of(
                &b.squares,
                cols,
                inner,
                "the transpose of the squares of the dosages of b",
            )?;
            let mut of_b = a_vector_of(0.0, num_values, &the_memory_for("Σy", num_values))?;
            let mut squares_of_b =
                a_vector_of(0.0, num_values, &the_memory_for("Σyy", num_values))?;
            sum_of(&a.called, &dosages_of_b, &mut of_b, "Σy")?;
            sum_of(
                &a.called,
                &squares_of_b_by_individual,
                &mut squares_of_b,
                "Σyy",
            )?;
            (of_b, squares_of_b)
        };
        Ok(TheSumsOfThePairs {
            num_individuals,
            products,
            of_a,
            of_b,
            squares_of_a,
            squares_of_b,
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
/// The six sums are whole numbers and each is held exactly in an `f64`,
/// so the only rounding is in the square and the division at the end.
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
        return f64::NAN;
    }
    above_the_line * above_the_line / (spread_of_a * spread_of_b)
}

/// The transpose of `matrix`, which holds `num_rows` rows of `num_cols`
/// values one after another and whose transpose holds `num_cols` rows of
/// `num_rows` values.
///
/// The three matrices of [`LdDosages`] hold exactly their variants times
/// their individuals, and so do the sums of a set of pairs: a matrix that
/// held more values than its rows times its columns would be transposed
/// up to that many, and one that held fewer would leave the rest of the
/// transpose at 0.
///
/// `what` names the transpose in the error of the memory.
///
/// # Errors
///
/// [`Error::LdNoMemory`] when this machine did not give the memory of the
/// transpose.
fn the_transpose_of(
    matrix: &[f64],
    num_rows: usize,
    num_cols: usize,
    what: &'static str,
) -> Result<Vec<f64>> {
    let values = matrix.len();
    let mut transposed = a_vector_of(0.0, values, &the_memory_for(what, values))?;
    if num_rows == 0 || num_cols == 0 {
        // A matrix with no row or no column has no value to transpose, and
        // neither of the two runs below is over a chunk of nothing.
        return Ok(transposed);
    }
    for (col, row_of_the_transpose) in transposed.chunks_exact_mut(num_rows).enumerate() {
        let values = row_of_the_transpose
            .iter_mut()
            .zip(matrix.chunks_exact(num_cols));
        for (value, row) in values {
            if let Some(found) = row.get(col) {
                *value = *found;
            }
        }
    }
    Ok(transposed)
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

/// The error of `values` values of one of the matrices of the r² that this
/// machine did not give the memory for, which `what` names.
fn the_memory_for(what: &'static str, values: usize) -> impl Fn() -> Error {
    move || Error::LdNoMemory { what, values }
}

/// How many values a matrix of `num_vars` variants of `num_individuals`
/// individuals holds, or `None` when they are more than the linear algebra
/// counts in.
fn the_values_of(num_vars: usize, num_individuals: usize) -> Option<usize> {
    num_vars
        .checked_mul(num_individuals)
        .filter(|values| *values <= MAX_VALUES_OF_THE_DOSAGES)
}

/// The `values` values of `matrix` from `from`, which is one matrix of a
/// run of variants of another.
fn the_values_from(matrix: &[f64], from: usize, values: usize) -> Vec<f64> {
    matrix.iter().skip(from).take(values).copied().collect()
}

/// The genotypes of the individuals of `individuals`, in the order they are
/// given, written into `chosen`.
///
/// `gts` is one row of the genotypes of a block, the alleles of one
/// individual after those of the individual before it, and `chosen` holds
/// one genotype for each index. [`LdDosages::of_block`] refuses an index
/// that is not an individual of the block before it reads a variant, so
/// every genotype is found; one that was not would leave that genotype of
/// `chosen` as it was.
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
        if let Some(theirs) = gts.get(from..to) {
            for (allele, theirs) in genotype.iter_mut().zip(theirs) {
                *allele = *theirs;
            }
        }
    }
}

/// The major allele frequency of the variant `counts` were counted for,
/// over its `called_alleles` called alleles, and `None` when it has none.
///
/// It is one division of the two counts as `f64`, as
/// `docs/specs/filters.md` has it, so that popnei gives one number for the
/// frequency wherever it is read.
fn the_major_allele_frequency(counts: &AlleleCounts, called_alleles: u32) -> Option<f64> {
    if called_alleles == 0 {
        return None;
    }
    let largest = counts.iter().copied().max().unwrap_or(0);
    Some(f64::from(largest) / f64::from(called_alleles))
}

/// The three values of each individual at one variant, and whether the
/// called genotypes of it hold two dosages at least.
///
/// `genotypes` holds one genotype of `of_a_genotype` alleles for each
/// individual, and `row`, `called` and `squares` one value each: the dosage
/// of the genotype, how many of its alleles are not `major`, with a 0 for a
/// genotype that has an allele missing; a 1 where the genotype was called
/// and a 0 where it was not; and the square of the dosage.
#[expect(
    clippy::arithmetic_side_effects,
    reason = "the dosage counts the alleles of one genotype, which are the ploidy, and \
              `LdDosages::of_block` refuses a ploidy above the 255 a u8 holds"
)]
fn the_dosages_of_a_variant(
    genotypes: &[i8],
    of_a_genotype: NonZeroUsize,
    major: i8,
    row: &mut [f64],
    called: &mut [f64],
    squares: &mut [f64],
) -> bool {
    let mut of_the_first_called = None;
    let mut has_variance = false;
    let values = row
        .iter_mut()
        .zip(called.iter_mut())
        .zip(squares.iter_mut())
        .zip(genotypes.chunks_exact(of_a_genotype.get()));
    for (((dosage_of, called_of), square_of), genotype) in values {
        let mut dosage = 0_u8;
        let mut missing = 0_u8;
        for allele in genotype {
            dosage += u8::from(*allele != major);
            missing |= u8::from(*allele == MISSING_ALLELE);
        }
        if missing != 0 {
            // A genotype with an allele missing has no dosage: it is a 0
            // in the dosages and in their squares, and a 0 in the called
            // genotypes takes it out of every sum of every pair its
            // variant is in.
            *dosage_of = 0.0;
            *called_of = 0.0;
            *square_of = 0.0;
            continue;
        }
        let value = f64::from(dosage);
        *dosage_of = value;
        *called_of = 1.0;
        *square_of = value * value;
        match of_the_first_called {
            None => of_the_first_called = Some(dosage),
            Some(first) => has_variance |= first != dosage,
        }
    }
    has_variance
}

/// The dosage that an entry of the matrix of the dosages holds.
#[expect(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "an entry of the matrix is a whole number from 0 to the ploidy, written by \
              `the_dosages_of_a_variant`, and `LdDosages::of_block` refuses a ploidy above \
              the 255 a u8 holds"
)]
fn the_dosage_of(value: f64) -> u8 {
    value as u8
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use super::{
        LdDosages, MAX_PLOIDY_OF_THE_DOSAGES, MAX_VALUES_OF_THE_DOSAGES, TheIndividualsThatDiffer,
        TheSumsOfThePairs, a_vector_of, r2_between, the_memory_for, the_values_of,
    };
    use crate::block::{Block, BlockReader};
    use crate::error::Error;
    use crate::io::vcf::{VcfOptions, VcfReader};
    use crate::variant::{MISSING_ALLELE, Needs};

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

    /// That two rows of one of the matrices hold the same values, within
    /// what the last bit of a whole number below 2^53 allows, which is 0.
    fn assert_the_values_are(found: &[f64], expected: &[f64], what: &str) {
        assert_eq!(
            found.len(),
            expected.len(),
            "{what}: the values are not as many"
        );
        for (at, (found, expected)) in found.iter().zip(expected).enumerate() {
            assert!(
                (found - expected).abs() < 1e-12,
                "{what}: the value {at} is {found} and not {expected}"
            );
        }
    }

    /// That two numbers are the same, within what the last bit of a
    /// division of two counts below 2^53 allows.
    fn assert_the_number_is(found: f64, expected: f64, what: &str) {
        assert!(
            (found - expected).abs() < 1e-15,
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
    fn a_matrix_this_machine_has_not_the_memory_for_is_an_error_and_not_the_end_of_the_process() {
        // The memory of every matrix of the r² is asked for with
        // `try_reserve_exact`, which gives it back as an error where
        // `vec![0.0; n]` would end the process, and which refuses a number
        // of values whose bytes this machine does not count before it asks
        // the allocator for anything.
        let values = usize::MAX;
        match a_vector_of(0.0_f64, values, &the_memory_for("the dosages", values)) {
            Err(Error::LdNoMemory {
                what,
                values: found,
            }) => {
                assert_eq!((what, found), ("the dosages", values));
            }
            Ok(given) => panic!("{} values of 8 bytes were given", given.len()),
            Err(other) => panic!("the memory failed with another error: {other:?}"),
        }
        let message = Error::LdNoMemory {
            what: "Σxy",
            values: 25,
        }
        .to_string();
        assert!(message.contains("Σxy"), "{message}");
        assert!(message.contains("25 values"), "{message}");
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
        let sums = TheSumsOfThePairs::of(&dosages, &dosages, 25).expect("the sums");
        for (of_a, of_b, expected, _) in THE_PAIRS_OF_THE_EXAMPLE {
            let found = [
                of_the_pair(&sums.num_individuals, 5, of_a, of_b),
                of_the_pair(&sums.of_a, 5, of_a, of_b),
                of_the_pair(&sums.of_b, 5, of_a, of_b),
                of_the_pair(&sums.products, 5, of_a, of_b),
                of_the_pair(&sums.squares_of_a, 5, of_a, of_b),
                of_the_pair(&sums.squares_of_b, 5, of_a, of_b),
            ];
            assert_the_sums_are(
                &found,
                &expected,
                &format!("the pair of the variants {of_a} and {of_b}"),
            );
            // The same pair the other way round has Σx and Σy, and Σxx and
            // Σyy, the other way round too.
            let [n, of_x, of_y, products, squares_of_x, squares_of_y] = expected;
            let back = [
                of_the_pair(&sums.num_individuals, 5, of_b, of_a),
                of_the_pair(&sums.of_a, 5, of_b, of_a),
                of_the_pair(&sums.of_b, 5, of_b, of_a),
                of_the_pair(&sums.products, 5, of_b, of_a),
                of_the_pair(&sums.squares_of_a, 5, of_b, of_a),
                of_the_pair(&sums.squares_of_b, 5, of_b, of_a),
            ];
            assert_the_sums_are(
                &back,
                &[n, of_y, of_x, products, squares_of_y, squares_of_x],
                &format!("the pair of the variants {of_b} and {of_a}"),
            );
        }
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
        let sums = TheSumsOfThePairs::of(&dosages, &dosages, 4).expect("the sums");
        assert_the_sums_are(
            &[
                of_the_pair(&sums.num_individuals, 2, 0, 1),
                of_the_pair(&sums.of_a, 2, 0, 1),
                of_the_pair(&sums.of_b, 2, 0, 1),
                of_the_pair(&sums.products, 2, 0, 1),
                of_the_pair(&sums.squares_of_a, 2, 0, 1),
                of_the_pair(&sums.squares_of_b, 2, 0, 1),
            ],
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
        let sums = TheSumsOfThePairs::of(&dosages, &dosages, 4).expect("the sums");
        assert_the_sums_are(
            &[
                of_the_pair(&sums.num_individuals, 2, 0, 1),
                of_the_pair(&sums.of_a, 2, 0, 1),
                of_the_pair(&sums.of_b, 2, 0, 1),
                of_the_pair(&sums.products, 2, 0, 1),
                of_the_pair(&sums.squares_of_a, 2, 0, 1),
                of_the_pair(&sums.squares_of_b, 2, 0, 1),
            ],
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

    /// One block with the 500 variants of the VCF at `path`, read as
    /// diploid and with the variants that failed their FILTER among them,
    /// which is what plink2 and pyNei were given.
    fn the_whole_of(path: &Path, needs: Needs) -> Block {
        let named = || path.display().to_string();
        let options = VcfOptions {
            ploidy: 2,
            only_passed: false,
            num_vars_per_block: Some(NUM_VARS_OF_A_REFERENCE),
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
            NUM_VARS_OF_A_REFERENCE,
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
        let block = the_whole_of(&the_reference_path("ld.vcf.gz"), Needs::GTS | Needs::ID);
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
        let block = the_whole_of(&many_vcf(), Needs::GTS);
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
}
