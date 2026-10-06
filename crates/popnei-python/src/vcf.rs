//! What a Python user reaches through `open_vcf`: a VCF with the options it
//! is read with, as a source of variants.
//!
//! [`VcfSource`] holds the path of a VCF and those options, and it reads the
//! header when it is built, so a file that is not a VCF fails at `open_vcf`.
//! Every pass over it opens the file again, which `source.rs` is what does:
//! the pass and the columns of its blocks are the same for a VCF and for a
//! vars file.

use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use pyo3::prelude::*;

use popnei::block::BlockReader;
use popnei::io::vcf::{
    VcfOptions, VcfReader, VcfWriteOptions, WriterSource, num_vars_per_block_of_write_vcf,
    ploidy_of_vcf,
};

use crate::errors::PyPopneiError;
use crate::source::{Blocks, OpenSource, PassCounts, blocks_of, count_of, source_of};
use crate::steps::Steps;
use crate::vars::write_the_pass;

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

    // Whether the variants that failed their FILTER are left out, which is
    // one of the two options a VCF is read with: the `repr` of a `Variants`
    // shows them, since two handles over one path that differ in them give
    // different variants.
    fn only_passed(&self) -> bool {
        self.options.only_passed
    }

    // The file the variants are read from, which the `repr` of a `Variants`
    // shows its user.
    fn path(&self) -> PathBuf {
        self.path.clone()
    }

    // One pass over the file, through the steps of `steps`: it is opened
    // again, and its blocks hold `fields` besides the genotypes,
    // `num_vars_per_block` variants each.
    #[pyo3(signature = (fields, num_vars_per_block, steps))]
    fn blocks(
        &self,
        py: Python<'_>,
        fields: Vec<String>,
        num_vars_per_block: Option<&Bound<'_, PyAny>>,
        steps: &Bound<'_, Steps>,
    ) -> Result<Blocks, PyPopneiError> {
        blocks_of(py, self, fields, num_vars_per_block, steps.get())
    }
}

impl OpenSource for VcfSource {
    fn path(&self) -> &Path {
        &self.path
    }

    fn ploidy(&self) -> usize {
        self.options.ploidy
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

    fn num_vars_per_block_of_the_vcf_writer(&self) -> Option<usize> {
        num_vars_per_block_of_write_vcf(WriterSource::Vcf {
            num_individuals: self.individuals.len(),
        })
    }
}

// The VCF at `path`, read with `ploidy` alleles in every genotype, or with
// the ploidy `ploidy_of_vcf` reads from the file when `ploidy` is `None`,
// and, when `only_passed` is true, without the variants that failed a
// filter. It reads the header, so the individuals are known when it
// returns.
#[pyfunction]
pub(crate) fn open_vcf(
    py: Python<'_>,
    path: PathBuf,
    ploidy: Option<&Bound<'_, PyAny>>,
    only_passed: bool,
) -> Result<VcfSource, PyPopneiError> {
    let ploidy = match ploidy {
        Some(ploidy) => count_of("ploidy", ploidy)?,
        None => py
            .detach(|| -> Result<_, popnei::Error> {
                let file = File::open(&path).map_err(|source| popnei::Error::FileNotOpened {
                    path: path.clone(),
                    source,
                })?;
                ploidy_of_vcf(BufReader::new(file))
            })
            // Every error of the search is of the file, a ploidy out of
            // range among them: it is the number of alleles of a genotype
            // the file holds, and not one the caller wrote.
            .map_err(|error| PyPopneiError::read_from_the_file(error, &path))?,
    };
    let options = VcfOptions {
        ploidy,
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

// Every variant of `source` into a VCF at `path`, through the steps of
// `steps`, bgzipped when the path ends in `.gz`, in any case, and plain
// otherwise.
// `source` is a VCF that `open_vcf` opened or a vars file that `open_vars`
// did, and what it gives back is the counts of the pass it made. The path
// is handled as `write_vars` handles its own: a path a file is at is
// refused, and the file of a call that failed is taken away. A `///`
// comment would become the `__doc__` of `popnei._core.write_vcf`, and what
// a Python user reads belongs to the package, which is the API.
#[pyfunction]
#[pyo3(signature = (source, path, steps))]
pub(crate) fn write_vcf(
    py: Python<'_>,
    source: &Bound<'_, PyAny>,
    path: PathBuf,
    steps: &Bound<'_, Steps>,
) -> Result<PassCounts, PyPopneiError> {
    let source = source_of(source)?;
    let steps = steps.get().of_a_pass()?;
    // `.gz` in any case of its two letters, `a.VCF.GZ` among them.
    let bytes = path.as_os_str().as_encoded_bytes();
    let options = VcfWriteOptions {
        bgzip: bytes
            .len()
            .checked_sub(3)
            .and_then(|start| bytes.get(start..))
            .is_some_and(|suffix| suffix.eq_ignore_ascii_case(b".gz")),
    };
    write_the_pass(
        py,
        source,
        &path,
        source.num_vars_per_block_of_the_vcf_writer(),
        &steps,
        |chain, sink| popnei::io::vcf::write_vcf(chain, sink, options),
    )
}
