"""The principal component analysis from Python, against pyNei.

`docs/specs/pca.md` has the two analyses, and both are here. `do_pca` takes
a pandas frame of individuals x traits. `do_pca_from_variants` takes a
`Variants` and makes one number of each variant for each individual, its
dosage, how many alleles of the genotype are not the major allele of the
variant. Each gives the projections of the individuals on the components,
the percentage of the variance each component holds and the weight of each
trait, or of each variant that was used, in each component.

The comparison both ask for is with pyNei at commit ef0ca6e, which
`pyproject.toml` names. For `do_pca` it is on `tests/reference/pca/iris.tsv`,
the 150 flowers x 4 measurements of pyNei's `test/datasets.py`, standardized
and not, and the four components of each are compared within 1e-9. For
`do_pca_from_variants` it is on `sim_missing.vcf`, the panel of 200
individuals and 1200 variants, and on `worked.vcf`, 5 individuals and 5
variants, the first 10 components of the one and the 3 of the other. A
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
from popnei import PCAResult, _core, do_pca, do_pca_from_variants, open_vcf
from pynei import vars_from_vcf
from pynei.pca import do_pca as pynei_do_pca
from pynei.pca import do_pca_from_variants as pynei_do_pca_from_variants
from pynei.var_filters import filter_by_maf

REFERENCE_PCA_DIR = Path(__file__).parent / "reference" / "pca"

# The two datasets of "How it is verified" of `docs/specs/pca.md` that both
# libraries are run on. The panel is pyNei's `test/gwas_reference/
# sim_missing.vars` written as a VCF: 200 individuals, 1200 variants of two
# alleles, three subpopulations, 7128 genotypes missing whole, every variant
# with variance. The worked example has 5 individuals and 5 variants, of
# which the ones at 0, 1 and 4 have variance; the variant at 4 has three
# alleles in the file and two among its called genotypes.
PANEL_VCF = REFERENCE_PCA_DIR / "sim_missing.vcf"
WORKED_VCF = REFERENCE_PCA_DIR / "worked.vcf"
# The worked example with a sixth variant of three alleles among its called
# genotypes, which needs `transform_to_biallelic`.
WORKED3_VCF = REFERENCE_PCA_DIR / "worked3.vcf"

PANEL_NUM_INDIVIDUALS = 200
PANEL_NUM_VARS = 1200
# What pyNei's `filter_by_maf` keeps of the panel at that threshold, from
# its `gather_filtering_stats` at ef0ca6e.
PANEL_MAF_THRESHOLD = 0.95
PANEL_VARS_KEPT_BY_THE_MAF_FILTER = 1175

# The tolerance of "How it is verified" of `docs/specs/pca.md`, read as an
# absolute difference. The largest number any comparison here holds is the
# 96.53 of a percentage of the iris run that centers nothing, whose
# projections reach 11.03; the standardized iris run reaches 72.96 and 3.31,
# and the one that centers without standardizing 92.46 and 3.80. Of the
# variants, the projections of the panel reach 17.10 and its percentages
# 7.61, and the worked example 3.03 and 76.74. At 96.53 a difference of 1e-9
# is the eleventh significant digit, and the reference files carry twelve.
# The two libraries differ by 1.0e-12 at most over the panel, where each
# takes another route to the same components.
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


def test_the_message_of_the_traits_with_no_variance_stops_at_ten():
    """Twelve traits with no variance: ten are named and two counted.

    A user of a table of hundreds of traits reads a message of one line and
    knows how many to take out, which is what pyNei's message does.
    """
    columns = {f"fixed{idx}": [1.0, 1.0, 1.0] for idx in range(12)}
    columns["a"] = [1.0, 2.0, 3.0]

    with pytest.raises(ValueError) as refusal:
        do_pca(pandas.DataFrame(columns))

    said = str(refusal.value)
    assert "12 of the 13 traits" in said
    assert all(f"`fixed{idx}`" in said for idx in range(10))
    assert "`fixed10`" not in said
    assert "`fixed11`" not in said
    assert "and 2 more" in said


@pytest.mark.parametrize(
    ("values", "what_it_says"),
    [
        ([1e308, 1.1e308, 1.2e308], "its mean is not finite"),
        ([1e154, -1e154, 3e154], "its standard deviation is not finite"),
        ([1e-200, 2e-200, 3e-200], "its standard deviation is 0"),
    ],
)
def test_a_trait_the_analysis_cannot_scale_is_refused_by_its_name(values, what_it_says):
    """Values too large or too small for the arithmetic of a float64.

    The values of a trait can sum above the largest float64, and then its
    mean is an infinity; the squares of its deviations can sum above it, and
    then its standard deviation is an infinity and the trait would become a
    column of zeros, looking like one with no variance; or those squares can
    all fall below the smallest float64 above 0, and then the deviation is 0
    although the values differ. The message names the trait, so that the
    user can scale that one or take it out.
    """
    data = pandas.DataFrame({"a": [1.0, 2.0, 3.0], "out_of_range": values})

    with pytest.raises(ValueError) as refusal:
        do_pca(data)

    said = str(refusal.value)
    assert "`out_of_range`" in said
    assert what_it_says in said


def test_a_table_where_no_trait_has_variance_is_refused():
    """There is no direction to give.

    Without standardizing a trait with no variance is no error, so a table
    of nothing but such traits reaches the decomposition, where every
    eigenvalue is 0. pyNei gives 0 for every projection and a percentage of
    NaN for every component.
    """
    data = pandas.DataFrame({"a": [1.0, 1.0, 1.0], "b": [2.0, 2.0, 2.0]})

    with pytest.raises(ValueError, match="no trait has variance"):
        do_pca(data, standardize_data=False)


def test_a_value_that_is_not_finite_is_refused_with_its_place():
    """An infinity, which pyNei lets through to numpy's decomposition.

    The message says which row and which trait it is at, because a table of
    hundreds of columns gives a user nowhere to look otherwise.
    """
    data = pandas.DataFrame({"a": [1.0, 2.0, 3.0], "b": [3.0, numpy.inf, 2.0]})

    with pytest.raises(ValueError, match="needs every value finite") as refusal:
        do_pca(data)

    assert "row 1, trait 1" in str(refusal.value)


def test_a_missing_value_of_a_nullable_dtype_is_refused_with_its_place():
    """pandas holds a missing value of a nullable dtype as its own NA.

    It reaches popnei as a NaN and is refused like any value that is not
    finite, with the row and the trait it is at, and not as an error about
    the dtype of the array.
    """
    data = pandas.DataFrame(
        {
            "a": pandas.array([1.0, 2.0, None], dtype="Float64"),
            "b": [3.0, 1.0, 2.0],
        }
    )

    with pytest.raises(ValueError, match="needs every value finite") as refusal:
        do_pca(data)

    assert "row 2, trait 0" in str(refusal.value)


def test_standardizing_a_table_that_is_not_centered_is_refused():
    """The deviation a trait is divided by is the one it has once centered.

    pyNei refuses the same pair.
    """
    data = pandas.DataFrame({"a": [1.0, 2.0, 3.0], "b": [3.0, 1.0, 2.0]})

    with pytest.raises(ValueError, match="standardized and not centered"):
        do_pca(data, center_data=False)


def test_a_table_of_fewer_than_two_rows_is_refused():
    """One individual has no variation for the components to hold.

    pyNei raises the error of the traits with no variance for it when it
    standardizes, and without standardizing divides by n - 1 = 0 and gives
    percentages of NaN.
    """
    data = pandas.DataFrame({"a": [1.0], "b": [3.0]})

    with pytest.raises(ValueError, match="2 rows at least") as refusal:
        do_pca(data)

    assert "the table is 1 x 2" in str(refusal.value)


def test_a_table_whose_products_are_not_finite_is_a_defect_of_popnei():
    """What the linear algebra refuses, which no argument of `do_pca` gives.

    The values of the table are finite and the products of the table with
    itself are not, which values of 1e200 give when nothing is centered or
    standardized to bring them down. Nothing a user writes is wrong here, so
    it is a ``RuntimeError`` and not a ``ValueError``, and it says which
    operation could not be done.
    """
    data = pandas.DataFrame({"a": [1e200, 2e200, 3e200], "b": [3e200, 1e200, 2e200]})

    with pytest.raises(RuntimeError, match="value that is not finite") as refusal:
        do_pca(data, center_data=False, standardize_data=False)

    assert "principal component analysis" in str(refusal.value)


def assert_the_first_components_are_pyneis(result, of_pynei, num_comps) -> None:
    """The first `num_comps` components of both libraries, within 1e-9.

    popnei gives the projections of every component that has variance and
    the weights of the first `num_prin_comps`, and pyNei gives both for the
    `min(individuals, variants)` components it has, so the comparison is of
    the first `num_comps` of each. The sign rule is applied to pyNei's
    numbers, each component's projections and weights together.
    """
    projections, princomps = the_signs_fixed(
        of_pynei.projections.to_numpy()[:, :num_comps],
        of_pynei.princomps.to_numpy()[:num_comps, :],
    )
    numpy.testing.assert_allclose(
        result.projections.to_numpy()[:, :num_comps],
        projections,
        rtol=0,
        atol=TOLERANCE,
    )
    numpy.testing.assert_allclose(
        result.explained_variance_percent.to_numpy()[:num_comps],
        of_pynei.explained_variance_percent.to_numpy()[:num_comps],
        rtol=0,
        atol=TOLERANCE,
    )
    numpy.testing.assert_allclose(
        result.princomps.to_numpy()[:num_comps, :], princomps, rtol=0, atol=TOLERANCE
    )


def test_the_first_ten_components_of_the_panel_are_pyneis():
    """The panel of 200 individuals and 1200 variants, with missing data.

    popnei gives 199 components, one fewer than pyNei's 200: centering takes
    one dimension out of the data and the last component has no variance.
    The weights are those of the first 10, the default, and every variant of
    the panel has variance, so all 1200 are the columns of `princomps`.
    """
    variants = open_vcf(PANEL_VCF)

    result = do_pca_from_variants(variants)

    of_pynei = pynei_do_pca_from_variants(vars_from_vcf(PANEL_VCF))
    assert result.projections.shape == (PANEL_NUM_INDIVIDUALS, 199)
    assert result.princomps.shape == (10, PANEL_NUM_VARS)
    assert_the_first_components_are_pyneis(result, of_pynei, 10)
    assert list(result.projections.index) == list(variants.individuals)
    assert list(result.projections.index) == list(of_pynei.projections.index)
    assert list(result.princomps.columns) == list(range(PANEL_NUM_VARS))
    # One component has one name in both frames: the width of the number is
    # that of the 199 components the projections have and not that of the 10
    # the weights are given for.
    assert list(result.princomps.index) == [f"PC00{number}" for number in range(10)]
    assert list(result.princomps.index) == list(result.projections.columns[:10])


def test_the_three_components_of_the_worked_example_are_pyneis():
    """5 individuals and 5 variants, of which 3 have variance.

    pyNei is given `transform_to_biallelic`, which it needs because it counts
    the alleles of the whole chunk and this file holds a variant with the
    alleles 0 and 1 and another with 0 and 2. It changes no dosage of a
    variant of two alleles, so the numbers are those of the run without it,
    which is what popnei does here.
    """
    variants = open_vcf(WORKED_VCF)

    result = do_pca_from_variants(variants, num_prin_comps=3)

    of_pynei = pynei_do_pca_from_variants(
        vars_from_vcf(WORKED_VCF), transform_to_biallelic=True
    )
    assert result.projections.shape == (5, 3)
    assert_the_first_components_are_pyneis(result, of_pynei, 3)
    assert list(result.princomps.columns) == [0, 1, 4]
    assert list(result.projections.columns) == ["PC0", "PC1", "PC2"]


def test_the_counts_of_the_pass_are_of_the_variants_the_steps_let_through():
    """`num_vars` counts the variants of the pass, used or not.

    The panel has no filter on it, so its count is every variant of the
    file, and `filtering` is empty.
    """
    result = do_pca_from_variants(open_vcf(PANEL_VCF), num_prin_comps=0)

    assert result.pass_stats.num_vars == PANEL_NUM_VARS
    assert result.pass_stats.filtering == {}


def test_a_filter_on_the_variants_is_counted_in_the_pass_stats():
    """The steps of the `Variants` run in both passes of the analysis.

    What the filter keeps is what pyNei's `filter_by_maf` keeps of the same
    file at the same threshold, and the variants of the pass are those. The
    weights are asked for, so the second pass runs through the steps as
    well: a second pass over other variants than the first is an error of
    the core, and one over the unfiltered variants would give other weights.
    """
    variants = open_vcf(PANEL_VCF)
    variants.filter_by_maf(PANEL_MAF_THRESHOLD)

    result = do_pca_from_variants(variants, num_prin_comps=10)

    of_pynei = pynei_do_pca_from_variants(
        filter_by_maf(vars_from_vcf(PANEL_VCF), max_allowed_maf=PANEL_MAF_THRESHOLD)
    )
    counts = result.pass_stats.filtering["maf"]
    assert counts.vars_processed == PANEL_NUM_VARS
    assert counts.vars_kept == PANEL_VARS_KEPT_BY_THE_MAF_FILTER
    assert result.pass_stats.num_vars == PANEL_VARS_KEPT_BY_THE_MAF_FILTER
    assert result.princomps.shape == (10, PANEL_VARS_KEPT_BY_THE_MAF_FILTER)
    assert list(result.princomps.columns) == list(of_pynei.princomps.columns)
    assert_the_first_components_are_pyneis(result, of_pynei, 10)


def test_no_weights_are_given_with_a_num_prin_comps_of_zero():
    """There is no second pass then, and the used variants are still there.

    `princomps` has no rows and its columns are the variants that were used,
    which is what says which variants the projections come from.
    """
    result = do_pca_from_variants(open_vcf(WORKED_VCF), num_prin_comps=0)

    assert result.princomps.shape == (0, 3)
    assert list(result.princomps.columns) == [0, 1, 4]
    assert result.projections.shape == (5, 3)


def test_more_components_than_there_are_gives_those_there_are():
    """The weights of 10 components of a table that has 3."""
    result = do_pca_from_variants(open_vcf(WORKED_VCF), num_prin_comps=10)

    assert result.princomps.shape == (3, 3)
    assert list(result.princomps.index) == list(result.projections.columns)


def test_the_positions_of_the_used_variants_are_numbers_with_a_sign():
    """The columns of the weights are whole numbers a user takes from.

    They are the positions of the variants that were used, and a user who
    asks what is one before each of them gets -1 for the variant at 0, as
    they do of pyNei's columns. Held as numbers with no sign, that first one
    would be 18446744073709551615.
    """
    result = do_pca_from_variants(open_vcf(WORKED_VCF), num_prin_comps=3)

    assert list(result.princomps.columns - 1) == [-1, 0, 3]


@pytest.mark.parametrize(
    ("given", "refusal", "said"),
    [
        (-1, ValueError, "-1"),
        (2**70, ValueError, "1180591620717411303424"),
        (2.5, TypeError, "2.5"),
        ("ten", TypeError, "'ten'"),
        (True, TypeError, "True"),
    ],
)
def test_a_num_prin_comps_that_counts_no_components_is_refused(given, refusal, said):
    """It says how many components the weights are given for.

    A whole number below 0 and one above what the machine counts are a
    ``ValueError``, and what is no whole number is a ``TypeError``: a truth
    value is a whole number in Python and would be the count 1 with nothing
    said. Each of them names the argument and what was given for it.
    """
    with pytest.raises(refusal) as refused:
        do_pca_from_variants(open_vcf(WORKED_VCF), num_prin_comps=given)

    assert "num_prin_comps" in str(refused.value)
    assert said in str(refused.value)


def test_do_pca_from_variants_refuses_what_is_not_a_variants():
    """A user who gives the path of the VCF, which is the easiest mistake.

    What it gave was the ``AttributeError`` of a `str` with no ``_source``.
    The refusal names the argument, says what was given and says where the
    variants come from, as the refusal of `write_vars` does.
    """
    with pytest.raises(TypeError, match="open_vcf") as refusal:
        do_pca_from_variants(str(WORKED_VCF))

    assert "variants" in str(refusal.value)
    assert repr(str(WORKED_VCF)) in str(refusal.value)


def test_a_variant_of_more_than_two_alleles_is_refused_by_its_position():
    """The dosage of a genotype has a meaning for a variant of two alleles.

    `worked3.vcf` is the worked example with a sixth variant whose called
    genotypes hold three alleles. The message says which variant it is and
    which argument reads it.
    """
    with pytest.raises(ValueError, match="transform_to_biallelic") as refusal:
        do_pca_from_variants(open_vcf(WORKED3_VCF))

    assert "the variant at the position 5" in str(refusal.value)


def test_transform_to_biallelic_reads_the_variant_of_three_alleles():
    """Every allele that is not the major one counts the same.

    popnei gives 3 components of `worked3.vcf` and pyNei 4, the last of
    which has a percentage below 1e-31 and is a component with no variance;
    the first 3 are the same in both.
    """
    variants = open_vcf(WORKED3_VCF)

    result = do_pca_from_variants(variants, transform_to_biallelic=True)

    of_pynei = pynei_do_pca_from_variants(
        vars_from_vcf(WORKED3_VCF), transform_to_biallelic=True
    )
    assert result.projections.shape == (5, 3)
    assert of_pynei.projections.shape == (5, 4)
    assert_the_first_components_are_pyneis(result, of_pynei, 3)
    assert list(result.princomps.columns) == [0, 1, 4, 5]
