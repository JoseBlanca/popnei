"""The configuration of the website of popnei's Python API.

Build it from the root of the repository, after `uv run maturin develop`,
with `uv run sphinx-build docs/api docs/api/_build`, and open
`docs/api/_build/index.html`.
"""

import popnei

project = "popnei"
author = "Jose Blanca"
copyright = "2026, Jose Blanca"
release = popnei.__version__
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
exclude_patterns = ["_build"]

# The docstrings write the name of an argument between single backticks,
# `vcf_path`, and a name of the package the same way, `Variants`. As a
# reference to a Python object the second is a link and the first, which
# no object is called, is set as code.
default_role = "py:obj"

# Every public name is imported from `popnei` itself, so the pages show
# `do_pca` and not `popnei.pca.do_pca`, which would read as a module a user
# has to import.
add_module_names = False
autodoc_member_order = "bysource"
autodoc_typehints = "signature"
autodoc_default_options = {"members": True, "show-inheritance": False}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
}

html_theme = "furo"
html_title = f"popnei {release}"
templates_path = ["_templates"]
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/scroll-start.html",
        "sidebar/home.html",
        "sidebar/navigation.html",
        "sidebar/scroll-end.html",
    ]
}


def _resolve_a_public_name(app, env, node, contnode):
    """The link of a name that the docstrings write as the package exports
    it, `popnei.Variants`, to the page of the object, which autodoc files
    under the module that defines it, `popnei.variant.Variants`.

    A name written without the package, `Variants.iter_blocks` in the
    docstring of another module, is resolved the same way.
    """
    if node.get("refdomain") != "py":
        return None
    target = node["reftarget"]
    head, _, rest = target.removeprefix("popnei.").partition(".")
    public = getattr(popnei, head, None)
    module = getattr(public, "__module__", None)
    if module is None or not module.startswith("popnei."):
        return None
    defined_at = f"{module}.{head}" + (f".{rest}" if rest else "")
    if defined_at == target:
        return None
    return env.domains["py"].resolve_xref(
        env, node["refdoc"], app.builder, node["reftype"], defined_at, node, contnode
    )


def setup(app):
    app.connect("missing-reference", _resolve_a_public_name)
