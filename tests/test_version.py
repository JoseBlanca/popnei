"""That the version a Python user sees is the one of the core crate.

It is the first test of the chain the Python side is made of: the core
crate, the binding crate `crates/popnei-python` that builds the module
`popnei._core`, and the package `python/popnei`. When it passes, the three
were built together and the package imports the module that was just
compiled.
"""

import importlib.metadata
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


def test_the_distribution_version_is_the_one_of_the_core_crate() -> None:
    """The version of the wheel, which is a second version of its own.

    maturin reads it from the manifest of the binding crate, and
    `popnei.__version__` comes the other way, through the compiled module
    from the core crate. The two are the same line of the workspace
    manifest only while both crates take `version.workspace = true`, so
    nothing but this catches a wheel published under another number.
    """
    assert importlib.metadata.version("popnei") == _version_of_the_core_crate()
