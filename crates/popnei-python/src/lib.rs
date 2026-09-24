//! The Python binding crate: it translates between Python and the core
//! crate `popnei` and holds no calculation of its own, as section 8 of
//! `docs/architecture.md` and `.claude/skills/coding/pyo3.md` say.
//!
//! What it builds is the private module `popnei._core`, which only the
//! Python package `python/popnei` imports. The package is what a user
//! sees: it carries the signatures of pyNei, the defaults and the
//! docstrings.

use pyo3::prelude::*;

mod dists;
mod errors;
mod gwas;
mod kinship;
mod ld;
mod pca;
mod pop_dists;
mod source;
mod stats;
mod steps;
mod vars;
mod vcf;

// The version of the core crate, `major.minor.patch`, which the Python
// package publishes as its own `__version__`. A `///` comment here would
// become the `__doc__` of `popnei._core.version`, and the documentation a
// Python user reads belongs to the package, which is the API.
#[pyfunction]
fn version() -> &'static str {
    popnei::version()
}

#[pymodule]
mod _core {
    // The defaults of `open_vcf` are the constants of the core crate, so
    // that the number a user gets when they say nothing is written in one
    // place. The Python package puts them in its signature.
    #[pymodule_export]
    const DEFAULT_PLOIDY: usize = popnei::io::vcf::DEFAULT_PLOIDY;
    #[pymodule_export]
    const DEFAULT_ONLY_PASSED: bool = popnei::io::vcf::DEFAULT_ONLY_PASSED;

    // The defaults of the statistics per population, which the Python
    // package puts in the signature of `calc_per_var_distribs`: how many
    // called genotypes a population needs at a variant to have a value
    // there, the histogram the values are counted in, and the major allele
    // frequency below which a variant is polymorphic.
    #[pymodule_export]
    const DEFAULT_MIN_NUM_INDIVIDUALS: u32 = popnei::stats::DEFAULT_MIN_NUM_INDIVIDUALS;
    #[pymodule_export]
    const DEFAULT_POLY_THRESHOLD: f64 = popnei::stats::DEFAULT_POLY_THRESHOLD;
    #[pymodule_export]
    const DEFAULT_HIST_RANGE: (f64, f64) = popnei::stats::DEFAULT_HIST_RANGE;
    #[pymodule_export]
    const DEFAULT_NUM_BINS: usize = popnei::stats::DEFAULT_NUM_BINS;
    #[pymodule_export]
    const DEFAULT_BIN_TYPE: &str = popnei::stats::DEFAULT_BIN_TYPE;

    // The two of `do_pca`, which are the constants of the core crate for
    // the same reason.
    #[pymodule_export]
    const DEFAULT_CENTER_DATA: bool = popnei::pca::DEFAULT_CENTER_DATA;
    #[pymodule_export]
    const DEFAULT_STANDARDIZE_DATA: bool = popnei::pca::DEFAULT_STANDARDIZE_DATA;

    // The two exceptions that carry the positions of the traits a
    // principal component analysis refused to `popnei.do_pca`, which names
    // those traits and raises the `ValueError` its user reads. A class made
    // with `create_exception!` is not a `#[pyclass]`, so it is added to the
    // module here.
    #[pymodule_init]
    fn add_the_exceptions(module: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
        use pyo3::prelude::PyModuleMethods as _;

        let py = module.py();
        module.add(
            "TraitsWithNoVariance",
            py.get_type::<super::errors::TraitsWithNoVariance>(),
        )?;
        module.add(
            "TraitOutOfRange",
            py.get_type::<super::errors::TraitOutOfRange>(),
        )?;
        // The pair of individuals of a kinship that have no variant called
        // in both of them, which `popnei.calc_kinship` catches to name the
        // two as the user knows them.
        module.add(
            "KinshipPairWithNoVariantCalled",
            py.get_type::<super::errors::KinshipPairWithNoVariantCalled>(),
        )
    }

    // The cap of `calc_rogers_huff_r2_matrix`, which is the constant of the
    // core: the matrix holds one r² for each pair of the variants of the
    // pass, so a pass of more than this many is refused instead of the
    // machine being asked for the memory of their square.
    #[pymodule_export]
    const DEFAULT_MAX_NUM_VARS: usize = popnei::ld::MAX_NUM_VARS_OF_THE_MATRIX;

    // The two of `do_pca_from_variants`, from the core as well.
    #[pymodule_export]
    const DEFAULT_TRANSFORM_TO_BIALLELIC: bool = popnei::pca::DEFAULT_TRANSFORM_TO_BIALLELIC;
    #[pymodule_export]
    const DEFAULT_NUM_PRIN_COMPS: usize = popnei::pca::DEFAULT_NUM_PRIN_COMPS;

    // Whether a mixed model of `calc_gwas` stands in for the denominator of
    // its test with the GRAMMAR-Gamma approximation when the user says
    // nothing. It is the exact denominator, since the approximation gives up
    // accuracy that grows with how strongly a panel is structured, and a
    // user asks for it.
    #[pymodule_export]
    const DEFAULT_USE_GRAMMAR_GAMMA_APPROX: bool = popnei::gwas::DEFAULT_USE_GRAMMAR_GAMMA_APPROX;

    #[pymodule_export]
    use super::dists::calc_pairwise_kosman_dists;
    #[pymodule_export]
    use super::gwas::calc_gwas;
    #[pymodule_export]
    use super::kinship::{calc_kinship, kinship_principal_components};
    #[pymodule_export]
    use super::ld::calc_rogers_huff_r2_matrix;
    #[pymodule_export]
    use super::pca::{pca, pca_of_variants};
    #[pymodule_export]
    use super::pop_dists::{calc_pop_dists, pop_dist_measures_that_have_a_value};
    #[pymodule_export]
    use super::source::Blocks;
    #[pymodule_export]
    use super::stats::{calc_per_individual_stats, calc_per_var_distribs};
    #[pymodule_export]
    use super::steps::Steps;
    #[pymodule_export]
    use super::vars::{VarsSource, open_vars, write_vars};
    #[pymodule_export]
    use super::vcf::{VcfSource, open_vcf};
    #[pymodule_export]
    use super::version;
}
