//! What a TypeScript user reaches through `doPca` and `doPcaFromVariants`:
//! the principal components of a table of numbers that the page holds,
//! individuals x traits, and those of the variants of a source.
//!
//! The table crosses as one `Float64Array`, row after row, which the code
//! wasm-bindgen generates copies into the memory of wasm before any code of
//! popnei runs, and the arrays of every result cross back the same way. The
//! variants are not a table that crosses: the core reads them block by block
//! from the readers this crate opens over the source and the steps of the
//! `Variants`, as it does for `write_vars`, and what comes back is of the
//! size of the individuals and of the variants that were used. Both analyses
//! are the core crate's and nothing of either is here.

use wasm_bindgen::prelude::wasm_bindgen;

use popnei::pca::{Pca, PcaOptions, VariantPcaOptions};

use crate::errors::JsPopneiError;
use crate::source::{OpenSource, PassCounts};
use crate::steps::{Steps, chain_of};

/// The principal components of a table, on their way to TypeScript.
///
/// Each array leaves the memory of wasm the first time it is asked for, and
/// the call after that gives nothing: the package reads each of them once,
/// into the object a user holds, and frees this. A copy left behind would
/// grow the memory of wasm, which never gives memory back, by the weights of
/// every column a second time.
#[wasm_bindgen]
pub struct PcaOfATable {
    num_comps: usize,
    projections: Option<Vec<f64>>,
    explained_variance_percent: Option<Vec<f64>>,
    num_prin_comps: usize,
    princomps: Option<Vec<f64>>,
}

#[wasm_bindgen]
impl PcaOfATable {
    /// How many components have variance, which is how many are given.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each row falls along each component, the rows of the table x
    /// `num_comps`, row after row.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }

    /// The variance of each component as a percentage of the variance of
    /// every component of the data, the ones with no variance counted, one
    /// number per component.
    pub fn explained_variance_percent(&mut self) -> Option<Vec<f64>> {
        self.explained_variance_percent.take()
    }

    /// How many components the weights are given for, which for a table is
    /// `num_comps`.
    #[must_use]
    pub fn num_prin_comps(&self) -> usize {
        self.num_prin_comps
    }

    /// The weight of each trait in each component, `num_prin_comps` x the
    /// traits of the table, row after row.
    pub fn princomps(&mut self) -> Option<Vec<f64>> {
        self.princomps.take()
    }
}

impl From<Pca> for PcaOfATable {
    /// The three arrays a table's analysis gives, without the two sides of
    /// the table, which the caller wrote, and without `used_cols`, which for
    /// a table is every trait.
    fn from(result: Pca) -> PcaOfATable {
        PcaOfATable {
            num_comps: result.num_comps,
            projections: Some(result.projections),
            explained_variance_percent: Some(result.explained_variance_percent),
            num_prin_comps: result.num_prin_comps,
            princomps: Some(result.princomps),
        }
    }
}

/// The principal components of the table `data`, `num_rows` rows of
/// `num_cols` values each, row after row, the rows being the individuals and
/// the columns the traits.
///
/// `center_data` takes the mean of each trait from it and `standardize_data`
/// then divides each trait by its standard deviation; standardizing without
/// centering is an error. No value may be missing.
///
/// The package checks that `data` holds exactly `num_rows` times `num_cols`
/// values before this is called, so the error the core has for a buffer that
/// is not of its table is not reached from TypeScript.
///
/// # Errors
///
/// When the analysis of a table of these two sides does not fit in the
/// memory of a page, when a value of the table is an infinity or a NaN, when
/// the table is to be standardized and not centered, when it has fewer than
/// 2 rows or no traits, when it is standardized and a trait has no variance,
/// when no trait of it has any, when the mean or the standard deviation of a
/// trait is not a number the analysis can use, and when the linear algebra
/// of the analysis could not be done.
#[wasm_bindgen]
pub fn pca(
    data: Vec<f64>,
    num_rows: usize,
    num_cols: usize,
    center_data: bool,
    standardize_data: bool,
) -> Result<PcaOfATable, JsPopneiError> {
    room_for_the_analysis_of_a_table(num_rows, num_cols)?;
    let options = PcaOptions {
        center: center_data,
        standardize: standardize_data,
    };
    Ok(popnei::pca::pca(&data, num_rows, num_cols, &options)?.into())
}

/// Whether the mean of each trait is taken from it when the caller says
/// nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_center_data() -> bool {
    popnei::pca::DEFAULT_CENTER_DATA
}

/// Whether each trait is divided by its standard deviation when the caller
/// says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_standardize_data() -> bool {
    popnei::pca::DEFAULT_STANDARDIZE_DATA
}

/// Whether every allele that is not the major one counts the same when the
/// caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_transform_to_biallelic() -> bool {
    popnei::pca::DEFAULT_TRANSFORM_TO_BIALLELIC
}

/// How many components the weights of the variants are given for when the
/// caller says nothing.
#[wasm_bindgen]
#[must_use]
pub fn default_num_prin_comps() -> usize {
    popnei::pca::DEFAULT_NUM_PRIN_COMPS
}

/// The principal components of the variants of a source, on their way to
/// TypeScript.
///
/// Each array leaves the memory of wasm the first time it is asked for, as
/// the arrays of a table's analysis do. The names of the individuals are not
/// here: the `Variants` a user gave read them from the header of the VCF or
/// the schema of the vars file when it was opened, and the package puts them
/// on the result.
#[wasm_bindgen]
pub struct PcaOfVariants {
    num_comps: usize,
    projections: Option<Vec<f64>>,
    explained_variance_percent: Option<Vec<f64>>,
    num_prin_comps: usize,
    princomps: Option<Vec<f64>>,
    used_vars: Option<Vec<u32>>,
    pass_stats: PassCounts,
}

#[wasm_bindgen]
impl PcaOfVariants {
    /// How many components have variance, which is how many are given.
    #[must_use]
    pub fn num_comps(&self) -> usize {
        self.num_comps
    }

    /// Where each individual falls along each component, the individuals of
    /// the source x `num_comps`, row after row.
    pub fn projections(&mut self) -> Option<Vec<f64>> {
        self.projections.take()
    }

    /// The variance of each component as a percentage of the variance of
    /// every component of the data, the ones with no variance counted, one
    /// number per component.
    pub fn explained_variance_percent(&mut self) -> Option<Vec<f64>> {
        self.explained_variance_percent.take()
    }

    /// How many components the weights are given for, which is
    /// `num_prin_comps` of the call, or the components there are when fewer
    /// were found, or 0 when none were asked for.
    #[must_use]
    pub fn num_prin_comps(&self) -> usize {
        self.num_prin_comps
    }

    /// The weight of each variant that was used in each component,
    /// `num_prin_comps` x the variants of `used_vars`, row after row.
    pub fn princomps(&mut self) -> Option<Vec<f64>> {
        self.princomps.take()
    }

    /// The position of each variant that was used, from 0, among the
    /// variants the pass gave, the ones with no variance included.
    pub fn used_vars(&mut self) -> Option<Vec<u32>> {
        self.used_vars.take()
    }

    /// How many variants the pass gave, used or not, and what each filter of
    /// it was given and kept. The two passes count the same and these are
    /// the counts of the first.
    #[must_use]
    pub fn pass_stats(&self) -> PassCounts {
        self.pass_stats.clone()
    }
}

/// The principal components of the variants of `source`, through the steps
/// of `steps`, with the weights of the first `num_prin_comps` components.
///
/// Both readers are opened here, over the same source and the same steps,
/// and lent to the core, which reads the variants of each in blocks: the
/// counts of the filters of the first are read when it returns, as
/// `docs/specs/filters.md` asks of every calculation. No reader is opened
/// for the second pass when no weights were asked for, and the core then
/// makes one pass.
///
/// # Errors
///
/// When the source cannot be read again, a wrong line of a VCF among the
/// causes; when a variant has more than two alleles among its called
/// genotypes and `transform_to_biallelic` is false; when the pass gives no
/// variant or no variant with variance; when a size of the dataset is beyond
/// what the analysis counts in; and when the linear algebra could not be
/// done.
pub(crate) fn pca_of_the_variants(
    source: &dyn OpenSource,
    num_individuals: usize,
    transform_to_biallelic: bool,
    num_prin_comps: usize,
    steps: Steps,
) -> Result<PcaOfVariants, JsPopneiError> {
    room_for_the_analysis_of_the_variants(num_individuals)?;
    let options = VariantPcaOptions {
        transform_to_biallelic,
        num_prin_comps,
    };
    // The source is asked for no size of block: the core puts a `reblock`
    // over each reader and chooses the size there, since the product of a
    // block is matrix work and a filter leaves blocks of uneven size.
    let mut first_pass = chain_of(source.reader(None)?, steps.steps())?;
    // The weights of a variant need the eigenvectors, which are known when
    // the first pass ends, so they come from a second pass over the same
    // variants. With none asked for there is no second reader and the source
    // is read once.
    let mut second_pass = if num_prin_comps > 0 {
        Some(chain_of(source.reader(None)?, steps.steps())?)
    } else {
        None
    };
    let result = popnei::pca::pca_of_variants(&mut first_pass, second_pass.as_mut(), &options)?;
    // The variants the pass gave, used or not, which is what the counts of a
    // pass say.
    let num_vars = u64::try_from(result.num_cols).map_err(|_| {
        JsPopneiError::Broken(format!(
            "the pass gave {num_vars} variants, more than the count of a pass holds",
            num_vars = result.num_cols
        ))
    })?;
    let pass_stats = PassCounts::of(num_vars, &first_pass.filtering_stats());
    Ok(PcaOfVariants {
        num_comps: result.num_comps,
        projections: Some(result.projections),
        explained_variance_percent: Some(result.explained_variance_percent),
        num_prin_comps: result.num_prin_comps,
        princomps: Some(result.princomps),
        used_vars: Some(the_positions_of_the_used_variants(&result.used_cols)?),
        pass_stats,
    })
}

/// The memory a wasm module addresses, 4 GiB, which is what a page holds of
/// everything popnei has open in it at once.
const MEMORY_OF_A_WASM_MODULE: u64 = 4 * 1024 * 1024 * 1024;

/// How much memory a principal component analysis holds at its peak, in
/// tenths of the square matrix it decomposes, which is 8 bytes per pair of
/// the smaller side of the data: that matrix, its eigenvectors, and the
/// workspace the eigendecomposition of faer allocates for itself.
///
/// The workspace is faer's own allocation and not one popnei asks for, and
/// an allocation that fails in wasm aborts, which is a trap that ends the
/// module where section 11 of `docs/architecture.md` asks for an `Error`. So
/// the size is measured and the analysis is refused before it starts.
///
/// Measured under node 24 on 22 September 2026, on a VCF of 2 variants read
/// with `numPrinComps` 1: 3000 individuals grew the memory of wasm to
/// 385220608 bytes, 5.35 times the 72000000 of their matrix; 9410
/// individuals ran, 9415 ended the module with `RuntimeError: unreachable`
/// after 173 ms, and 9415 individuals have a matrix of 709137800 bytes, of
/// which 4 GiB is 6.05. So the analysis holds about 6 times its matrix, and
/// popnei counts 6.1 of them, which takes 9381 of a side at most, 29 below
/// the smallest number that trapped. A table is decomposed by the same code
/// of the same library over the same matrix, so the same count holds for the
/// smaller of its two sides.
const TENTHS_OF_THE_MATRIX_THE_ANALYSIS_HOLDS: u64 = 61;

/// Which side of the data the matrix that is decomposed is of, for the
/// message of an analysis that does not fit.
#[derive(Clone, Copy)]
enum TheSquareOf {
    /// The individuals of a dataset of variants.
    Individuals,
    /// The rows of a table that has fewer rows than traits.
    Rows,
    /// The traits of a table that has fewer traits than rows.
    Traits,
}

impl TheSquareOf {
    /// What the side is called in the message.
    fn named(self) -> &'static str {
        match self {
            Self::Individuals => "individuals",
            Self::Rows => "rows",
            Self::Traits => "traits",
        }
    }
}

/// That the memory of wasm takes the principal components of a dataset whose
/// smaller side is `side` of `what`.
///
/// What the analysis holds is
/// [`TENTHS_OF_THE_MATRIX_THE_ANALYSIS_HOLDS`] tenths of the square matrix of
/// that side, and a page holds 4 GiB of everything at once. What this does
/// not know is what the tab already holds, the bytes of the file among them,
/// so a page with little left can still run out; what it stops is the
/// dataset that cannot fit however empty the tab is, which is the one that
/// ends the module with no message.
///
/// # Errors
///
/// When the analysis of that many does not fit in the memory of a page.
fn room_for_the_square_of(side: usize, what: TheSquareOf) -> Result<(), JsPopneiError> {
    // A count that is beyond what these multiplications hold is a dataset
    // that is far beyond the memory of a page, so every one of them
    // saturates instead of being checked: what the number then says is the
    // largest the arithmetic holds, and the analysis is refused either way.
    let pairs = u64::try_from(side).unwrap_or(u64::MAX);
    let wanted = pairs
        .saturating_mul(pairs)
        .saturating_mul(8)
        .saturating_mul(TENTHS_OF_THE_MATRIX_THE_ANALYSIS_HOLDS)
        / 10;
    if wanted <= MEMORY_OF_A_WASM_MODULE {
        return Ok(());
    }
    Err(JsPopneiError::NoMemory(format!(
        "the principal components of {side} {named} hold about {gigabytes} GB, the \
         {named} x {named} matrix of the analysis, its eigenvectors and the workspace \
         of the eigendecomposition, and a page holds at most 4 GB of everything at a \
         time. Data this large is analysed by a program outside the browser, popnei in \
         Python among them.",
        named = what.named(),
        gigabytes = wanted.div_ceil(1000 * 1000 * 1000)
    )))
}

/// That the memory of wasm takes the principal components of the variants of
/// `num_individuals` individuals, which are decomposed over the individuals
/// x individuals matrix.
///
/// # Errors
///
/// When the analysis of that many individuals does not fit in the memory of
/// a page.
#[wasm_bindgen]
pub fn room_for_the_analysis_of_the_variants(num_individuals: usize) -> Result<(), JsPopneiError> {
    room_for_the_square_of(num_individuals, TheSquareOf::Individuals)
}

/// That the memory of wasm takes the principal components of a table of
/// `num_rows` rows and `num_cols` traits.
///
/// The matrix that is decomposed is the square of the smaller of the two
/// sides, so a table of 150 rows and 4 traits is decomposed over 4 x 4 and
/// one of 4 rows and 150 traits over 4 x 4 as well.
///
/// # Errors
///
/// When the analysis of that table does not fit in the memory of a page.
#[wasm_bindgen]
pub fn room_for_the_analysis_of_a_table(
    num_rows: usize,
    num_cols: usize,
) -> Result<(), JsPopneiError> {
    if num_rows <= num_cols {
        room_for_the_square_of(num_rows, TheSquareOf::Rows)
    } else {
        room_for_the_square_of(num_cols, TheSquareOf::Traits)
    }
}

/// The positions of the variants that were used as the `Uint32Array` they
/// cross in.
///
/// # Errors
///
/// When a position is above 4294967295, which no pass of wasm reaches: a
/// `usize` is 32 bits there, and the core refuses a pass of more variants
/// than one counts.
fn the_positions_of_the_used_variants(used_cols: &[usize]) -> Result<Vec<u32>, JsPopneiError> {
    used_cols
        .iter()
        .map(|position| {
            u32::try_from(*position).map_err(|_| {
                JsPopneiError::NotInJavaScript(format!(
                    "the variant at the position {position} is above 4294967295, the \
                     largest position the array of the variants that were used holds"
                ))
            })
        })
        .collect()
}
