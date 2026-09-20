//! popnei: population genetics over the variants of a dataset, the
//! successor of the Python library pyNei.
//!
//! This is the core crate, where every calculation lives. It is pure Rust:
//! no Python and no JavaScript reach it, and it builds for the two targets
//! of WebAssembly that popnei ships, `wasm32-unknown-unknown` for a web
//! application and `wasm32-unknown-emscripten` for the wheel of pyodide.
//! Python and TypeScript call it through a binding crate each, as section 8
//! of `docs/architecture.md` lays out.
//!
//! A variant is one site of the genome with its alleles and the genotype of
//! every individual at it. The `variant` module has it, with the trait that
//! anything giving variants implements; the modules that read them, hold
//! them in blocks and calculate over them are being written, and
//! `docs/architecture.md` has their order.

#![forbid(unsafe_code)]

pub mod error;
pub mod variant;

pub use error::{Error, Result};

/// The version of this crate, `major.minor.patch`, as its manifest gives
/// it.
///
/// The Python and the TypeScript packages publish it as their own version,
/// so that a user who reports a result names the code that gave it.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::version;
    use flate2::Compression;
    use flate2::read::GzDecoder;
    use flate2::write::GzEncoder;
    use std::io::{Read, Write};

    /// The literal is the version of the workspace manifest, and it changes
    /// with it.
    #[test]
    fn version_is_the_one_of_the_manifest() {
        assert_eq!(version(), "0.1.0");
    }

    /// flate2 is here for the VCF reader, which decompresses gzip with it.
    /// The test reads back what it wrote through the default backend of
    /// flate2, miniz_oxide, the one that builds for both wasm targets.
    #[test]
    fn gzip_is_read_back_as_it_was_written() {
        let header = b"##fileformat=VCFv4.3\n";
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(header).unwrap();
        let gzipped = encoder.finish().unwrap();

        let mut read_back = Vec::new();
        GzDecoder::new(gzipped.as_slice())
            .read_to_end(&mut read_back)
            .unwrap();
        assert_eq!(read_back, header);
    }
}
