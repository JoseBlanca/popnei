//! The Python binding crate: it translates between Python and the core
//! crate `popnei` and holds no calculation of its own, as section 8 of
//! `docs/architecture.md` and `.claude/skills/coding/pyo3.md` say.
//!
//! What it builds is the private module `popnei._core`, which only the
//! Python package `python/popnei` imports. The package is what a user
//! sees: it carries the signatures of pyNei, the defaults and the
//! docstrings.

use pyo3::prelude::*;

mod errors;
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

    #[pymodule_export]
    use super::source::Blocks;
    #[pymodule_export]
    use super::stats::calc_per_var_distribs;
    #[pymodule_export]
    use super::steps::Steps;
    #[pymodule_export]
    use super::vars::{VarsSource, open_vars, write_vars};
    #[pymodule_export]
    use super::vcf::{VcfSource, open_vcf};
    #[pymodule_export]
    use super::version;
}
