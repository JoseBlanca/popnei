//! What a TypeScript user reaches through `calcKinship`: how much more of
//! their genome every pair of the individuals of a source share than two
//! individuals drawn at random from the same panel would.
//!
//! The calculation is the core's, `popnei::kinship::calc_kinship`, and what
//! this module does is the translation that section 11 of
//! `docs/architecture.md` leaves to a binding crate. It builds the chain of
//! readers of the pass from the steps of the `Variants`, turns the names of
//! the individuals a user asked the kinship of into their positions among
//! those the pass gives, and keeps that chain while the calculation runs so
//! that the counts of its filters can be read when it returns.
//!
//! [`KinshipOfVariants`] is the result on its way out. It lives in the
//! memory of wasm, which the garbage collector of JavaScript does not see,
//! so the package frees it as soon as its parts are read. The matrix is
//! moved out of it as it is read and not cloned: it is 8 bytes for each pair
//! of individuals, 800 MB at 10000 of them, and the memory of wasm never
//! gives back what it grew by.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::block::{Block, BlockReader};
use popnei::filters::{FilteringStats, resolve_individuals};
use popnei::kinship::calc_kinship;
use popnei::variant::{ChromTable, Needs};

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts};
use crate::steps::{Steps, chain_of};

/// The kinship of every pair of individuals, the names of those individuals
/// and the counts of the pass it was taken over.
///
/// The matrix is individuals x individuals, row after row, and it is
/// symmetric: the entry of the individuals `i` and `j` is the value at
/// `i * N + j`, and the value at `j * N + i` is the same one.
#[wasm_bindgen]
pub struct KinshipOfVariants {
    /// How many individuals the matrix has on each of its two sides.
    num_individuals: usize,
    /// How many variants had variance among these individuals and were
    /// used, which is what the matrix was built from.
    num_vars: usize,
    /// The matrix, and `None` once it was given to JavaScript: it leaves
    /// the memory of wasm as it is read, so the copy that crosses is the
    /// only one.
    matrix: Option<Vec<f64>>,
    /// The names of the individuals the matrix is of, in its order, and
    /// `None` once they were given to JavaScript.
    individuals: Option<Vec<String>>,
    /// The counts of the pass, which the package turns into the `passStats`
    /// of the result.
    counts: PassCounts,
}

#[wasm_bindgen]
impl KinshipOfVariants {
    /// How many individuals the matrix has on each of its two sides.
    #[must_use]
    pub fn num_individuals(&self) -> usize {
        self.num_individuals
    }

    /// How many variants the matrix was built from: those that had variance
    /// among these individuals. A variant whose called genotypes all have
    /// one dosage is in no sum and in no denominator, and the counts of the
    /// pass say how many variants the steps gave, used or not.
    #[must_use]
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// The matrix, individuals x individuals row after row, or `undefined`
    /// when it was read already.
    pub fn matrix(&mut self) -> Option<Vec<f64>> {
        self.matrix.take()
    }

    /// The names of the individuals, in the order of the rows and the
    /// columns of the matrix, or `undefined` when they were read already.
    pub fn individuals(&mut self) -> Option<Vec<String>> {
        self.individuals.take()
    }

    /// How many variants the pass gave, used or not, and what each filter of
    /// it was given and kept.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.counts.clone()
    }
}

/// The kinship of the individuals of `source`, over the variants the steps
/// of `steps` keep.
///
/// `individuals` are the names of the ones the matrix is of, in the order it
/// has them, and `None` is every individual the pass gives, in its order.
/// Every frequency, mean and denominator is of those individuals, so the
/// kinship of some of them is not the rows and the columns of the kinship of
/// the whole panel. `transform_to_biallelic` says that a variant of more
/// than two alleles among its called genotypes is read with every allele
/// that is not the major one counting the same, which is what
/// `pca_of_variants` also takes.
///
/// The chain of readers of the pass stays here, lent to the core, so that
/// the counts of its filters can be read when the calculation returns. How
/// many variants that chain gave is counted here as well, by
/// [`TheVariantsCounted`]: the core's `Kinship` carries the variants that
/// were used and not the ones the pass gave, and the counts of a pass are of
/// what the steps let through.
///
/// The source is asked for no size of block: the core puts a `reblock` over
/// the reader and chooses the size there, since the product of a block is
/// matrix work and a filter leaves blocks of uneven size.
///
/// # Errors
///
/// When a name of `individuals` is of nobody the pass gives, is there twice,
/// or the list is empty; when the source cannot be read, a wrong line of a
/// VCF among the causes; when a variant has more than two alleles among its
/// called genotypes and `transform_to_biallelic` is false; when the pass
/// gives no variant or no variant with variance; when two individuals have
/// no variant called in both; when a size of the dataset is beyond what the
/// calculation counts in; and when the linear algebra could not be done.
pub(crate) fn kinship_of_the_variants(
    source: &dyn OpenSource,
    individuals: Option<Vec<String>>,
    transform_to_biallelic: bool,
    steps: Steps,
) -> Result<KinshipOfVariants, JsPopneiError> {
    let chain = chain_of(source.reader(None)?, steps.steps())?;
    let mut counted = TheVariantsCounted::over(chain);
    // The names the pass gives, which are the source's own when no step is
    // a filter of individuals and the kept ones in the order they were
    // named when one is. They are read before the calculation borrows the
    // chain, so the matrix and the names cannot be of two different passes.
    let of_the_pass = counted.individuals().to_vec();
    // A name that is of nobody is refused before the source is read: the
    // rule and its message are the core's, the one a filter of individuals
    // is given its names by.
    let positions = match individuals.as_deref() {
        Some(names) => Some(resolve_individuals(names, &of_the_pass)?),
        None => None,
    };
    let kinship = calc_kinship(&mut counted, positions.as_deref(), transform_to_biallelic)?;
    let counts = PassCounts::of(counted.num_vars(), &counted.filtering_stats());
    // The names of the matrix are the ones that were asked for, in the
    // order they were asked in, which is the order the core gives the rows
    // in; with no name at all they are every individual of the pass.
    let of_the_matrix = individuals.unwrap_or(of_the_pass);
    Ok(KinshipOfVariants {
        num_individuals: kinship.num_individuals,
        num_vars: kinship.num_vars,
        matrix: Some(kinship.matrix),
        individuals: Some(of_the_matrix),
        counts,
    })
}

/// The chain of a pass with a count of the variants it gives, which is what
/// the counts of a pass say and what the kinship of the core does not carry.
///
/// `popnei::kinship::Kinship` has `num_vars`, the variants that had variance
/// and were used, and a variant with none is in no sum and in no
/// denominator, so it is not the number the `passStats` of the result holds.
/// Every other calculation of the core gives that number away with its
/// result, the principal components of the variants as `num_cols`, and the
/// kinship gives no block of its pass to this crate for it to be counted
/// anywhere else.
///
/// It reads no block and changes none: every method is the reader's below.
struct TheVariantsCounted<R: BlockReader> {
    reader: R,
    num_vars: u64,
}

impl<R: BlockReader> TheVariantsCounted<R> {
    /// The reader with its count at 0.
    fn over(reader: R) -> TheVariantsCounted<R> {
        TheVariantsCounted {
            reader,
            num_vars: 0,
        }
    }

    /// How many variants the reader has given so far.
    fn num_vars(&self) -> u64 {
        self.num_vars
    }
}

impl<R: BlockReader> BlockReader for TheVariantsCounted<R> {
    fn next_block(&mut self) -> popnei::Result<Option<Block>> {
        let block = self.reader.next_block()?;
        if let Some(block) = block.as_ref() {
            // A pass of wasm reads a file that is in the memory of the tab,
            // which addresses 2^32 bytes, so the count is nowhere near what
            // a `u64` holds; it saturates rather than wrap, because a count
            // that went round would be a smaller number than the truth and
            // nothing would say so.
            self.num_vars = self
                .num_vars
                .saturating_add(block.num_vars.try_into().unwrap_or(u64::MAX));
        }
        Ok(block)
    }

    fn individuals(&self) -> &[String] {
        self.reader.individuals()
    }

    fn ploidy(&self) -> usize {
        self.reader.ploidy()
    }

    fn chroms(&self) -> &ChromTable {
        self.reader.chroms()
    }

    fn set_needs(&mut self, needs: Needs) {
        self.reader.set_needs(needs);
    }

    fn filtering_stats(&self) -> Vec<(&'static str, FilteringStats)> {
        self.reader.filtering_stats()
    }
}
