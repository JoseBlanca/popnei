"""That the version a Python user sees is the one of the core crate.

It is the first test of the chain the Python side is made of: the core
crate, the binding crate `crates/popnei-python` that builds the module
`popnei._core`, and the package `python/popnei`. When it passes, the three
were built together and the package imports the module that was just
compiled.
"""

import tomllib
from pathlib import Path

import popnei

REPO_DIR = Path(__file__).parent.parent


def _version_of_the_core_crate() -> str:
    """Read the version of `crates/popnei` from the workspace manifest.

    The crate takes `version.workspace = true`, so the version is written
    once, in `[workspace.package]` of the `Cargo.toml` of the repository.
    """
    manifest = tomllib.loads((REPO_DIR / "Cargo.toml").read_text())
    return manifest["workspace"]["package"]["version"]


def test_the_package_version_is_the_version_of_the_core_crate() -> None:
    assert popnei.__version__ == _version_of_the_core_crate()
