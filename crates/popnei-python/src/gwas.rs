//! The association study of the variants of a source on its way to Python:
//! one pass over them, and one row for each with the effect of the variant
//! on the trait, how uncertain that effect is and its p-value.
//!
//! The calculation is the core's, [`popnei::gwas::calc_gwas`], and what this
//! module does is what `kinship.rs` does for the kinship: it builds the
//! chain of readers of the pass from the steps of the `Variants`, lends it
//! to the core with the interpreter released, and reads the counts of the
//! filters from that chain when the call is over, since no block of the pass
//! reaches this crate.
//!
//! Which individuals are tested, and their phenotype and their design, are
//! the package's: "Which individuals are tested, and the design" of
//! `docs/specs/gwas.md` puts the refusals about a name and about a frame in
//! the layer that holds the names, and the core is given the positions of
//! the tested individuals with the phenotype and the design already built.
//! So what crosses here is three arrays that are read together row by row,
//! and the core refuses an order that is not the source's, a phenotype that
//! is not a number and a design whose columns are not independent.
//!
//! The three arrays are copied into memory of Rust before the pass starts,
//! and not read where numpy holds them: the interpreter is released for the
//! whole pass, and what goes into that closure has to be owned by it. The
//! copy is one position and one trait for each tested individual and one row
//! of the design for each, which for a panel of 10000 individuals with three
//! covariates is 80 KB, 80 KB and 320 KB, against the variants x individuals
//! the pass itself reads.
//!
//! A mixed model brings a fourth, the kinship of those individuals, and that
//! one is lent where numpy holds it and not copied: it grows with the square
//! of the panel, 800 MB for 10000 individuals, and the array is the caller's
//! and out of reach of Python for as long as the pass runs, which is what
//! `pca.rs` relies on for its table.

use numpy::{
    IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods as _,
};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use popnei::block::BlockReader;
use popnei::gwas::{Gwas, GwasInput, TestType, TraitType};

use crate::errors::{PyPopneiError, raise_a_ctrl_c_before_numpy_is_called};
use crate::source::{ChromColumn, OpenSource, PassCounts, chrom_column, id_column, source_of};
use crate::steps::{Step, Steps, chain_of};

/// The null model on its way to Python: which of the four models was
/// fitted, which test was made of every variant, the effect of each column
/// of the design, what the model left unexplained, the variance of the
/// random effect, the heritability and how many individuals were tested.
///
/// The names of the columns of the design are the package's, which has them
/// from the covariates its user gave.
type NullModelForPython<'py> = (
    &'static str,
    &'static str,
    Bound<'py, PyArray1<f64>>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
    usize,
);

/// The rows of the result on their way to Python, column by column: the
/// chromosome, the position and the id of each variant, each of them `None`
/// when the source has no such column, and then the four numbers of every
/// variant, which hold NaN for one that has no answer.
type StatsForPython<'py> = (
    Option<Bound<'py, PyTuple>>,
    Option<Bound<'py, PyArray1<u64>>>,
    Option<Bound<'py, PyTuple>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<f64>>,
);

/// What a study gives Python: the null model, the rows of the variants,
/// whether the GRAMMAR-Gamma approximation was used and the counts of the
/// pass.
type GwasForPython<'py> = (
    NullModelForPython<'py>,
    StatsForPython<'py>,
    bool,
    PassCounts,
);

/// What each filter of a chain was given and kept, under the kind of the
/// filter and in the order of the chain, which the package turns around for
/// its user.
type FilteringCounts = Vec<(&'static str, u64, u64)>;

// The study of the variants of `source` that the steps of `steps` keep, over
// the individuals at the positions `individuals` among the ones the pass
// gives, with `phenotype` their trait and `design` their design, one row of
// its columns each and the column of ones of the intercept first. A `///`
// comment here would become the `__doc__` of `popnei._core.calc_gwas`, and
// what a Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, individuals, phenotype, design, trait_name, test_name, kinship, use_grammar_gamma_approx, transform_to_biallelic, steps))]
#[expect(
    clippy::too_many_arguments,
    reason = "what a study is given in `docs/specs/gwas.md`: the source and its steps, \
              the three arrays of the tested individuals that are read together row by \
              row, the kinship of those individuals, and the three things a user asked \
              for. A struct of them would be built in Python, element by element"
)]
pub(crate) fn calc_gwas<'py>(
    py: Python<'py>,
    source: &Bound<'_, PyAny>,
    individuals: PyReadonlyArray1<'py, u64>,
    phenotype: PyReadonlyArray1<'py, f64>,
    design: PyReadonlyArray2<'py, f64>,
    trait_name: &str,
    test_name: Option<String>,
    kinship: Option<PyReadonlyArray2<'py, f64>>,
    use_grammar_gamma_approx: bool,
    transform_to_biallelic: bool,
    steps: &Bound<'_, Steps>,
) -> Result<GwasForPython<'py>, PyPopneiError> {
    let source = source_of(source)?;
    // The two names of a trait and the two of a test are the core's, which
    // is where a name that is of neither is refused: one list of them
    // serves both packages, and the message a user reads is the same in
    // each. Which tests the model of the study has is the core's too, and
    // it is `the_model_and_the_test` that answers it.
    let trait_type = TraitType::of_name(trait_name)?;
    let test = test_name.as_deref().map(TestType::of_name).transpose()?;
    let tested = the_positions(&individuals)?;
    let trait_values = the_values("phenotype", &phenotype)?;
    // The columns of the design are its second dimension, so the core is
    // given a count that the array itself says and not one this crate was
    // told beside it.
    let (_, num_coefs) = design.as_array().dim();
    // The layout is asked of the array itself and not of `as_slice`, which
    // takes an array that lies column after column as well: the core reads
    // the design row after row, one row for each tested individual, and a
    // design that lies the other way would be another matrix of the same
    // numbers.
    if !design.is_c_contiguous() {
        return Err(PyPopneiError::ArrayNotContiguous { name: "design" });
    }
    let design_values = design
        .as_slice()
        .map_err(|_| PyPopneiError::ArrayNotContiguous { name: "design" })?
        .to_vec();
    // The kinship crosses as the matrix of the tested individuals alone,
    // cut and ordered by the package, which is where the names of the
    // individuals are: the core reads it row after row beside the design
    // and holds no name to cut it by. It is lent to the pass where numpy
    // holds it and not copied, as `pca.rs` lends its table: the array is
    // owned by the caller of this function and nothing of Python can reach
    // it while the interpreter is released, and a copy would be 800 MB for
    // 10000 individuals where the three arrays beside it are kilobytes.
    let kinship_values = the_kinship(kinship.as_ref())?;
    let steps = steps.get().of_a_pass()?;
    // A Ctrl-C that was pending when this was called is raised here, before
    // the file is opened.
    py.check_signals()?;
    let path = source.path().to_path_buf();
    let input = GwasInput {
        phenotype: &trait_values,
        trait_type,
        design: &design_values,
        num_coefs,
        kinship: kinship_values,
        test,
        // The approximation is refused by the core, which says two
        // different things about it, that a study with no kinship has no
        // denominator to approximate and that popnei has not written the
        // approximation of the one a mixed model has. Both packages give it
        // as the user wrote it, so both messages are the same in each.
        use_grammar_gamma_approx,
        individuals: &tested,
        transform_to_biallelic,
    };
    // The whole source is read inside this one call, and the test of a block
    // is matrix work, so the interpreter is released for all of it: the
    // threads of rayon deadlock on a caller that holds it. A Ctrl-C that
    // arrives meanwhile is raised when the call is over and not between two
    // blocks, as it is in `Blocks::__next__`: the loop over the blocks is
    // the core's, and a study that is interrupted loses only itself, since
    // it writes no file and the `Variants` is as it was.
    let calculated = py.detach(|| over_the_source(source, &steps, &input));
    // The file of the source goes into every error of the pass, and
    // `errors.rs` is what leaves it out of the message of the ones that are
    // of an argument a user wrote.
    let (result, filtering) =
        calculated.map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    // The Ctrl-C that arrived while the interpreter was released is raised
    // before numpy is called: the first array of a process imports the C API
    // of numpy, that import fails with the exception that is pending, and
    // the numpy crate panics when it does, which a user cannot catch.
    raise_a_ctrl_c_before_numpy_is_called(py)?;
    // The result is taken apart rather than read field by field, so that
    // every column of it is given to numpy without a copy: the four of them
    // are 32 bytes for each variant, 32 MB for a million.
    let Gwas {
        num_vars,
        null_model,
        allele_freq,
        beta,
        se,
        p_value,
        used_grammar_gamma_approx,
        chroms,
        chrom_table,
        poss,
        ids,
    } = result;
    // The variants of the study go to Python as the `u64` every count of
    // popnei is there: a `usize` is 32 bits in WebAssembly and 64 natively,
    // and what a user reads does not depend on that.
    let num_vars = u64::try_from(num_vars).map_err(|_| {
        PyPopneiError::broken_of_the_file(
            format!("the study read {num_vars} variants, more than the count of a pass holds"),
            &path,
        )
    })?;
    let chroms = match chroms.as_deref() {
        Some(numbers) => {
            let column = ChromColumn::of(numbers, &chrom_table, &path)?;
            Some(chrom_column(py, &column, &path)?)
        }
        None => None,
    };
    let stats = (
        chroms,
        poss.map(|poss| poss.into_pyarray(py)),
        ids.as_deref().map(|ids| id_column(py, ids)).transpose()?,
        allele_freq.into_pyarray(py),
        beta.into_pyarray(py),
        se.into_pyarray(py),
        p_value.into_pyarray(py),
    );
    let null_model = (
        null_model.model.name(),
        null_model.test.name(),
        null_model.covariate_effects.into_pyarray(py),
        null_model.residual_variance,
        null_model.genetic_variance,
        null_model.heritability,
        null_model.num_individuals,
    );
    Ok((
        null_model,
        stats,
        used_grammar_gamma_approx,
        (num_vars, filtering),
    ))
}

/// One pass over `source` through `steps`, and the study of the variants it
/// gives.
///
/// The interpreter is already released here, so nothing of Python is
/// touched: what comes back is Rust.
///
/// # Errors
///
/// When the source cannot be opened or read, when the pass gives no variant,
/// and when the core refuses the study: a model that is not written yet, a
/// test the model has not, individuals that are not in the order the source
/// has them, a phenotype or a value of the design that is not a finite
/// number, a design whose columns are not independent or that leaves nothing
/// to measure the uncertainty of a variant from, a variant of more than two
/// alleles while `transform_to_biallelic` is false, and an operation of the
/// linear algebra that did not run.
fn over_the_source(
    source: &dyn OpenSource,
    steps: &[Step],
    input: &GwasInput<'_>,
) -> popnei::Result<(Gwas, FilteringCounts)> {
    // The source is opened at the size of its own blocks, since the core
    // puts a `reblock` over whatever it is given.
    let reader = source.reader(None)?;
    // The chain of the pass stays here, lent to the core, so that the counts
    // of its filters can be read when the call is over: the loop over the
    // blocks is the core's, and no block of it reaches this crate.
    let mut chain = chain_of(reader, steps)?;
    let result = popnei::gwas::calc_gwas(&mut chain, None::<&mut Box<dyn BlockReader>>, input)?;
    let filtering = chain
        .filtering_stats()
        .into_iter()
        .map(|(kind, stats)| (kind, stats.vars_processed, stats.vars_kept))
        .collect();
    Ok((result, filtering))
}

/// The positions of the tested individuals as the core counts them.
///
/// They cross as an array of whole numbers, which is what one value per
/// individual crosses as, and the core takes them as the whole numbers this
/// build indexes memory with, 32 bits in WebAssembly.
///
/// # Errors
///
/// [`PyPopneiError::ArrayNotContiguous`] when the array does not lie value
/// after value, and [`PyPopneiError::Broken`] when a position is more than
/// this build indexes with, which the package cannot reach: it builds the
/// positions by walking the individuals the pass gives.
fn the_positions(individuals: &PyReadonlyArray1<'_, u64>) -> Result<Vec<usize>, PyPopneiError> {
    individuals
        .as_slice()
        .map_err(|_| PyPopneiError::ArrayNotContiguous {
            name: "individuals",
        })?
        .iter()
        .map(|position| {
            usize::try_from(*position).map_err(|_| PyPopneiError::Broken {
                message: format!(
                    "the individual at the position {position} was asked to be tested, and \
                     this build counts the individuals of a source to {largest}",
                    largest = usize::MAX
                ),
                path: None,
            })
        })
        .collect()
}

/// The kinship of the tested individuals as the core reads it, row after
/// row, and `None` for a study with no random effect.
///
/// Whether it holds one row and one column for each tested individual, and
/// whether every value of it is a finite number, are the core's to refuse:
/// it says which cell a value that is not finite is in, and the package
/// checks neither, since `Kinship.__post_init__` refuses a matrix that is
/// not square and one that holds what is no number.
///
/// # Errors
///
/// [`PyPopneiError::ArrayNotContiguous`] when the array does not lie row
/// after row, which the package makes it before the call. A matrix that
/// lies column after column would be read as its own transpose, which for
/// a kinship is the same matrix; the layout is asked for all the same,
/// because `as_slice` takes a strided view of another matrix altogether.
fn the_kinship<'a>(
    kinship: Option<&'a PyReadonlyArray2<'_, f64>>,
) -> Result<Option<&'a [f64]>, PyPopneiError> {
    match kinship {
        None => Ok(None),
        Some(values) => {
            if !values.is_c_contiguous() {
                return Err(PyPopneiError::ArrayNotContiguous { name: "kinship" });
            }
            Ok(Some(values.as_slice().map_err(|_| {
                PyPopneiError::ArrayNotContiguous { name: "kinship" }
            })?))
        }
    }
}

/// The values of an array of one number per tested individual, copied into
/// memory of Rust so that the pass can own them while the interpreter is
/// released.
///
/// # Errors
///
/// [`PyPopneiError::ArrayNotContiguous`] when the array does not lie value
/// after value, which the package makes it before the call.
fn the_values(
    name: &'static str,
    values: &PyReadonlyArray1<'_, f64>,
) -> Result<Vec<f64>, PyPopneiError> {
    Ok(values
        .as_slice()
        .map_err(|_| PyPopneiError::ArrayNotContiguous { name })?
        .to_vec())
}
