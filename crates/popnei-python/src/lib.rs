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
mod pca;
mod source;
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
        )
    }

    // The two of `do_pca_from_variants`, from the core as well.
    #[pymodule_export]
    const DEFAULT_TRANSFORM_TO_BIALLELIC: bool = popnei::pca::DEFAULT_TRANSFORM_TO_BIALLELIC;
    #[pymodule_export]
    const DEFAULT_NUM_PRIN_COMPS: usize = popnei::pca::DEFAULT_NUM_PRIN_COMPS;

    #[pymodule_export]
    use super::dists::calc_pairwise_kosman_dists;
    #[pymodule_export]
    use super::pca::{pca, pca_of_variants};
    #[pymodule_export]
    use super::source::Blocks;
    #[pymodule_export]
    use super::steps::Steps;
    #[pymodule_export]
    use super::vars::{VarsSource, open_vars, write_vars};
    #[pymodule_export]
    use super::vcf::{VcfSource, open_vcf};
    #[pymodule_export]
    use super::version;
}
