"""The principal component analysis from Python, against pyNei.

`docs/specs/pca.md` has the two analyses. What is here is the one of a
table, `do_pca`, which takes a pandas frame of individuals x traits and
gives the projections of the individuals on the components, the percentage
of the variance each component holds and the weight of each trait in each
component.

The comparison it asks for is with pyNei at commit ef0ca6e, which
`pyproject.toml` names: both libraries run on `tests/reference/pca/iris.tsv`,
the 150 flowers x 4 measurements of pyNei's `test/datasets.py`, standardized
and not, and the four components of each are compared within 1e-9. A
component multiplied by -1 is the same component, and pyNei gives whichever
its decomposition gave, so the sign rule of "What both analyses compute" of
the spec is applied to pyNei's numbers here before they are compared: in
each component the projection of the largest absolute value is made
positive, and the weights of that component take the same sign.

The run without centering is compared with pyNei alone: "How it is verified"
of the spec says that no program outside the project gives a number for it,
since R's `prcomp` with `center=FALSE` divides by n - 1 a sum of squares
that was not centered.
"""

from pathlib import Path

import numpy
import pandas
import pytest
from popnei import PCAResult, _core, do_pca
from pynei.pca import do_pca as pynei_do_pca

REFERENCE_PCA_DIR = Path(__file__).parent / "reference" / "pca"

# The tolerance of "How it is verified" of `docs/specs/pca.md`. The
# projections of iris reach 3.3, so it is read as an absolute difference.
TOLERANCE = 1e-9

# The table of `test_pca_refuses_traits_with_no_variance` of pyNei, whose
# trait `fixed` is the same number in every row.
NO_VARIANCE_TABLE = pandas.DataFrame(
    {"a": [1.0, 2.0, 3.0], "fixed": [5.0, 5.0, 5.0], "b": [3.0, 1.0, 2.0]}
)


@pytest.fixture
def iris() -> pandas.DataFrame:
    """The 150 flowers x 4 measurements both libraries are run on.

    The first column of the file is the index, the number of the flower, and
    the other four are named after the measurement each one holds.
    """
    return pandas.read_csv(REFERENCE_PCA_DIR / "iris.tsv", sep="\t", index_col=0)


def the_signs_fixed(
    projections: numpy.ndarray, princomps: numpy.ndarray
) -> tuple[numpy.ndarray, numpy.ndarray]:
    """pyNei's numbers under the sign rule of `docs/specs/pca.md`.

    In each component the projection of the largest absolute value is made
    positive, the first of them when two individuals have the same absolute
    value, and the weights of that component are multiplied by the same -1
    or 1. It is what popnei does with its own numbers, so that the result
    does not depend on which library did the decomposition.
    """
    projections = numpy.array(projections)
    princomps = numpy.array(princomps)
    for component in range(projections.shape[1]):
        largest = numpy.argmax(numpy.abs(projections[:, component]))
        if projections[largest, component] < 0:
            projections[:, component] *= -1
            princomps[component, :] *= -1
    return projections, princomps


def assert_it_is_pyneis(result: PCAResult, of_pynei) -> None:
    """The three tables of both libraries, within 1e-9 after the sign rule.

    Every component pyNei gives is compared: on the tables here popnei gives
    as many as pyNei does, so a component that popnei left out for having no
    variance would fail the first assertion.
    """
    projections, princomps = the_signs_fixed(
        of_pynei.projections.to_numpy(), of_pynei.princomps.to_numpy()
    )
    assert result.projections.shape == of_pynei.projections.shape
    numpy.testing.assert_allclose(
        result.projections.to_numpy(), projections, rtol=0, atol=TOLERANCE
    )
    numpy.testing.assert_allclose(
        result.explained_variance_percent.to_numpy(),
        of_pynei.explained_variance_percent.to_numpy(),
        rtol=0,
        atol=TOLERANCE,
    )
    numpy.testing.assert_allclose(
        result.princomps.to_numpy(), princomps, rtol=0, atol=TOLERANCE
    )


@pytest.mark.parametrize("standardize_data", [True, False])
def test_the_four_components_of_iris_are_pyneis(iris, standardize_data):
    """The default run and the one that centers without standardizing."""
    result = do_pca(iris, standardize_data=standardize_data)
    of_pynei = pynei_do_pca(iris, standarize_data=standardize_data)
    assert_it_is_pyneis(result, of_pynei)


def test_iris_without_centering_is_pyneis(iris):
    """The run that centers nothing, which only pyNei gives a number for."""
    result = do_pca(iris, center_data=False, standardize_data=False)
    of_pynei = pynei_do_pca(iris, center_data=False, standarize_data=False)
    assert_it_is_pyneis(result, of_pynei)


def test_the_result_carries_the_names_of_the_frame(iris):
    """The index and the columns of the table are those of the result.

    The components are named `PC0` to `PC3`, and the result of a table has
    no pass over a source of variants behind it, so its `pass_stats` is
    `None`.
    """
    result = do_pca(iris)
    names = ["PC0", "PC1", "PC2", "PC3"]
    assert list(result.projections.index) == list(iris.index)
    assert list(result.projections.columns) == names
    assert list(result.explained_variance_percent.index) == names
    assert list(result.princomps.index) == names
    assert list(result.princomps.columns) == list(iris.columns)
    assert result.pass_stats is None


def test_the_names_of_ten_components_and_more_have_a_zero_on_the_left():
    """A table of 12 traits, whose components are `PC00` to `PC11`.

    The width of the number is that of the count of the components, as
    pyNei's `_create_pc_names` gives it, so that the names of a hundred
    components sort as their numbers do.
    """
    numbers = numpy.random.default_rng(0).normal(size=(15, 12))
    data = pandas.DataFrame(numbers, columns=[f"trait{idx}" for idx in range(12)])
    result = do_pca(data)
    assert list(result.projections.columns) == [f"PC{idx:02d}" for idx in range(12)]
    assert_it_is_pyneis(result, pynei_do_pca(data))


def test_the_private_module_names_an_array_that_is_not_contiguous(iris):
    """What a user who calls `popnei._core` themselves reads.

    The package makes every array lie in memory row after row before the
    call, so nothing a user writes reaches this message; when it is read, it
    names the argument and what makes an array the core can take.
    """
    by_columns = numpy.asfortranarray(iris.to_numpy())

    with pytest.raises(ValueError, match="ascontiguousarray") as refusal:
        _core.pca(by_columns, True, True)

    assert "`data`" in str(refusal.value)


def test_a_trait_with_no_variance_is_refused_by_its_name():
    """Standardizing divides by a standard deviation of 0.

    The user is told how many traits have none and which they are, by the
    name each one has in the frame, so that they can take them out.
    """
    with pytest.raises(ValueError, match="no variance"):
        do_pca(NO_VARIANCE_TABLE)
    with pytest.raises(ValueError, match="fixed"):
        do_pca(NO_VARIANCE_TABLE)


def test_a_trait_with_no_variance_is_no_error_without_standardizing():
    """It gets a weight of 0, and the table has 2 components and not 3.

    Centering takes one of the three dimensions of the rows out, so the
    third component has no variance and is not given, which is where popnei
    differs from pyNei: pyNei gives 3 components, the last of them numerical
    noise. The numbers are in "How it is verified" of `docs/specs/pca.md`.
    """
    result = do_pca(NO_VARIANCE_TABLE, standardize_data=False)
    assert result.projections.shape == (3, 2)
    numpy.testing.assert_allclose(
        result.explained_variance_percent.to_numpy(), [75.0, 25.0], rtol=0, atol=1e-9
    )
    numpy.testing.assert_allclose(
        result.princomps["fixed"].to_numpy(), [0.0, 0.0], rtol=0, atol=1e-9
    )
