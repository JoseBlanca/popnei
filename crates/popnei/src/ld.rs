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

use std::num::NonZeroUsize;

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
    /// How many individuals each variant has a value for, which are the
    /// ones [`LdDosages::of_block`] was given.
    num_individuals: usize,
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
    /// order they are given, and an empty slice is every individual of it.
    /// The major allele of each variant is that of
    /// [`the_major_allele`](crate::variant::the_major_allele) over those
    /// individuals alone, so the dosages of a population are counted from
    /// the allele that population was called most often at.
    ///
    /// # Errors
    ///
    /// A block that does not pass [`Block::check`], one with variants and
    /// no genotypes, an index that is not an individual of the block, a
    /// block whose genotypes hold more than
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
        for individual in individuals {
            if *individual >= block.num_individuals {
                return Err(Error::LdIndividualNotInTheDataset {
                    individual: *individual,
                    num_individuals: block.num_individuals,
                });
            }
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
            num_individuals,
            dosages: vec![0.0; values],
            called: vec![0.0; values],
            squares: vec![0.0; values],
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
        self.num_individuals
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
        let from = first.saturating_mul(self.num_individuals);
        let values = num_vars.saturating_mul(self.num_individuals);
        Ok(LdDosages {
            num_vars,
            num_individuals: self.num_individuals,
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
        let from = var.checked_mul(self.num_individuals)?;
        let to = from.checked_add(self.num_individuals)?;
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
    use super::{LdDosages, MAX_PLOIDY_OF_THE_DOSAGES, MAX_VALUES_OF_THE_DOSAGES, the_values_of};
    use crate::block::Block;
    use crate::error::Error;
    use crate::variant::MISSING_ALLELE;

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
    fn a_block_of_no_variant_gives_dosages_of_no_variant() {
        let block = block_of(&[], 6, 2);
        let dosages = LdDosages::of_block(&block, &[]).expect("the dosages");
        assert_eq!((dosages.num_vars(), dosages.num_individuals()), (0, 6));
        assert!(dosages.dosages(0).is_none(), "a variant has dosages");
        assert_eq!(dosages.maf(0), None);
    }
}
