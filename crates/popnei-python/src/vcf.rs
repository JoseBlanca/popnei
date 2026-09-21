//! What a Python user reaches through `open_vcf`: a VCF with the options it
//! is read with, as a source of variants.
//!
//! [`VcfSource`] holds the path of a VCF and those options, and it reads the
//! header when it is built, so a file that is not a VCF fails at `open_vcf`.
//! Every pass over it opens the file again, which `source.rs` is what does:
//! the pass and the columns of its blocks are the same for a VCF and for a
//! vars file.

use std::path::{Path, PathBuf};

use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::io::vcf::{VcfOptions, VcfReader};

use crate::errors::PyPopneiError;
use crate::source::{Blocks, OpenSource, blocks_of, count_of};

// A VCF that was opened: its path, the options it is read with, and the
// individuals its header named. A `///` here would become the `__doc__` of
// the class, and what a Python user reads belongs to the package, which is
// the API.
#[pyclass(frozen, module = "popnei._core")]
pub(crate) struct VcfSource {
    path: PathBuf,
    options: VcfOptions,
    individuals: Vec<String>,
}

#[pymethods]
impl VcfSource {
    // The names of the individuals, in the order of the columns of the
    // VCF.
    fn individuals(&self) -> Vec<String> {
        self.individuals.clone()
    }

    // How many alleles the genotype of one individual holds.
    fn ploidy(&self) -> usize {
        self.options.ploidy
    }

    // One pass over the file: it is opened again, and its blocks hold
    // `fields` besides the genotypes, `num_vars_per_block` variants each.
    #[pyo3(signature = (fields, num_vars_per_block))]
    fn blocks(
        &self,
        py: Python<'_>,
        fields: Vec<String>,
        num_vars_per_block: Option<&Bound<'_, PyAny>>,
    ) -> Result<Blocks, PyPopneiError> {
        blocks_of(py, self, fields, num_vars_per_block)
    }
}

impl OpenSource for VcfSource {
    fn path(&self) -> &Path {
        &self.path
    }

    fn reader(
        &self,
        num_vars_per_block: Option<usize>,
    ) -> Result<Box<dyn BlockReader>, popnei::Error> {
        // The VCF reader cuts its blocks where it is asked to, so the size
        // the caller wants is one of the options it is built with.
        let options = VcfOptions {
            num_vars_per_block,
            ..self.options
        };
        Ok(Box::new(VcfReader::from_path(&self.path, options)?))
    }
}

// The VCF at `path`, read with `ploidy` alleles in every genotype and, when
// `only_passed` is true, without the variants that failed a filter. It
// reads the header, so the individuals are known when it returns.
#[pyfunction]
pub(crate) fn open_vcf(
    py: Python<'_>,
    path: PathBuf,
    ploidy: &Bound<'_, PyAny>,
    only_passed: bool,
) -> Result<VcfSource, PyPopneiError> {
    let options = VcfOptions {
        ploidy: count_of("ploidy", ploidy)?,
        only_passed,
        num_vars_per_block: None,
    };
    let individuals = py
        .detach(|| -> Result<_, popnei::Error> {
            // The header is read when the reader is built and no variant
            // is. Nothing here asks for a block, so a file whose blocks
            // would need more memory than this machine gives is opened all
            // the same and its individuals read, and the size of its blocks
            // is the user's to choose at `iter_blocks`.
            let reader = VcfReader::from_path(&path, options)?;
            Ok(reader.individuals().to_vec())
        })
        .map_err(|error| PyPopneiError::of_the_file(error, &path))?;
    Ok(VcfSource {
        path,
        options,
        individuals,
    })
}
