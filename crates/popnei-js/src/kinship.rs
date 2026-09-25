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

use popnei::block::BlockReader;
use popnei::filters::resolve_individuals;
use popnei::kinship::{calc_kinship, principal_components_of};

use crate::errors::JsPopneiError;
use crate::source::{Consumer, OpenSource, PassCounts, the_run_of};
use crate::steps::{Steps, chain_of};

/// The kinship of every pair of individuals, the names of those individuals
/// and the counts of the pass it was taken over.
///
/// The matrix is individuals x individuals, row after row, and it is
/// symmetric: the entry of the individuals `i` and `j` is the value at
/// `i * N + j`, and the value at `j * N + i` is the same one.
#[wasm_bindgen]
pub struct KinshipOfVariants {
    /// How many variants had variance among these individuals and were
    /// used, which is what the matrix was built from. It is an `f64`, the
    /// number of JavaScript, as every count `PassCounts` carries is: the
    /// core holds it in a `u64` and a count above 2^53 is more variants
    /// than any source has.
    num_vars: f64,
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
    /// How many variants the matrix was built from: those that had variance
    /// among these individuals. A variant whose called genotypes all have
    /// one dosage is in no sum and in no denominator, and the counts of the
    /// pass say how many variants the steps gave, used or not.
    #[must_use]
    pub fn num_vars(&self) -> f64 {
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
/// many variants it gave, used or not, is `Kinship::num_vars_given`, which
/// is the `num_vars` of those counts; `Kinship::num_vars` is the variants
/// that had variance and were used, and a variant with none is in no sum and
/// in no denominator.
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
    the_run_of(source, &Consumer::Kinship, |run| {
        let mut chain = chain_of(source.reader(run, None)?, steps.steps())?;
        // The names the pass gives, which are the source's own when no step is
        // a filter of individuals and the kept ones in the order they were
        // named when one is. They are read before the calculation borrows the
        // chain, so the matrix and the names cannot be of two different
        // passes.
        let of_the_pass = chain.individuals().to_vec();
        // A name that is of nobody is refused before the source is read: the
        // rule and its message are the core's, the one a filter of individuals
        // is given its names by.
        let positions = match individuals.as_deref() {
            Some(names) => Some(resolve_individuals(names, &of_the_pass)?),
            None => None,
        };
        // The names of the matrix are the ones that were asked for, in the
        // order they were asked in, which is the order the core gives the rows
        // in; with no name at all they are every individual of the pass. They
        // are built before the call because the pair that has no variant
        // called in both is named with them.
        let of_the_matrix = individuals.unwrap_or(of_the_pass);
        let kinship = calc_kinship(&mut chain, positions.as_deref(), transform_to_biallelic)
            .map_err(|error| under_the_names_of_the_individuals(error, &of_the_matrix))?;
        let counts = PassCounts::of(kinship.num_vars_given, &chain.filtering_stats());
        Ok(KinshipOfVariants {
            num_vars: kinship.num_vars as f64,
            matrix: Some(kinship.matrix),
            individuals: Some(of_the_matrix),
            counts,
        })
    })
}

/// Where each individual of a kinship falls along the directions in which
/// the panel varies most, on their way to TypeScript.
///
/// The projections leave the memory of wasm the first time they are asked
/// for and the call after that gives nothing, as the matrix of
/// [`KinshipOfVariants`] does: the package reads them once, into the object
/// a user holds, and frees this.
#[wasm_bindgen]
pub struct PcsOfAKinship {
    /// How many components were given, which is the `num_pcs` that were
    /// asked for or the components the matrix has when it has fewer.
    num_comps: usize,
    /// The individuals x `num_comps` projections, row after row, and `None`
    /// once they were given to JavaScript.
    projections: Option<Vec<f64>>,
}

#[wasm_bindgen]
impl PcsOfAKinship {
    /// How many components were given: a kinship measures a pair against
    /// the average pair of the panel, which takes one direction out of it
    /// when no genotype is missing, and both reference panels of
    /// `docs/specs/kinship.md` have 199 components and not 200.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each individual falls along each component, the individuals x
    /// `num_comps` row after row, or `undefined` when they were read
    /// already.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }
}

/// The principal components of the kinship `matrix`, the individuals of
/// `num_individuals` x the same individuals row after row, `num_pcs` of them
/// at most.
///
/// With `lambda_j` the eigenvalues of the matrix from the largest and `u_j`
/// its eigenvectors, the component `j` is `u_j * sqrt(lambda_j)`: where each
/// individual falls along the direction in which the panel varies the `j`th
/// most. A component whose eigenvalue is not above the tolerance of
/// `docs/specs/pca.md` is not given, and in every component the projection
/// of the largest absolute value is positive, which is the rule that makes
/// the two backends of the eigendecomposition and the three builds of
/// popnei give one answer.
///
/// The matrix is the one a `Kinship` of the package holds, a calculated one
/// or one a user built, so the package is what has checked that it is
/// square, symmetric and finite before this is called.
///
/// # Errors
///
/// When the matrix holds no individual, and when the eigendecomposition
/// could not be done.
#[wasm_bindgen]
pub fn kinship_principal_components(
    matrix: Vec<f64>,
    num_individuals: usize,
    num_pcs: usize,
) -> Result<PcsOfAKinship, JsPopneiError> {
    let num_values = matrix.len();
    if num_individuals.checked_mul(num_individuals) != Some(num_values) {
        // The `Kinship` of the package holds one value for each pair of its
        // individuals, which its constructor is what checks, so a caller
        // that arrives here wrote the matrix itself. Python says the same
        // of the same matrix, with a `ValueError`.
        return Err(JsPopneiError::Refused(format!(
            "the matrix of a kinship of {num_individuals} individuals holds \
             one value for each pair of them, and this one holds \
             {num_values}"
        )));
    }
    // The matrix is taken over and the eigendecomposition writes the
    // eigenvectors over it: a `Kinship` built here would carry two counts
    // nobody gave and the matrix would be copied to protect a kinship that
    // is thrown away.
    let pcs = principal_components_of(matrix, num_individuals, num_pcs)?;
    Ok(PcsOfAKinship {
        num_comps: pcs.num_comps,
        projections: Some(pcs.projections),
    })
}

/// `error` with the two individuals that have no variant called in both of
/// them under their names, and every other error as it is.
///
/// The core names the positions the two have among the individuals of the
/// kinship, which is what it has: a user drops a name from the panel, and
/// with `individuals` on the call those positions are not even the ones the
/// file has. A position the names do not reach is left as the core wrote it
/// rather than named wrongly; the core takes both from the matrix it built
/// of these names, so none is.
fn under_the_names_of_the_individuals(
    error: popnei::Error,
    individuals: &[String],
) -> JsPopneiError {
    let popnei::Error::KinshipPairWithNoVariantCalled {
        one,
        other,
        num_vars_of_one,
        num_vars_of_other,
    } = error
    else {
        return JsPopneiError::Core(error);
    };
    match (individuals.get(one), individuals.get(other)) {
        (Some(one), Some(other)) => JsPopneiError::PairWithNoVariantCalled {
            one: one.clone(),
            other: other.clone(),
            num_vars_of_one,
            num_vars_of_other,
        },
        _ => JsPopneiError::Core(popnei::Error::KinshipPairWithNoVariantCalled {
            one,
            other,
            num_vars_of_one,
            num_vars_of_other,
        }),
    }
}
