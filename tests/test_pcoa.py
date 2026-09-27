"""The principal coordinates of distances from Python, against pyNei.

`docs/specs/pca.md` has the analysis, in "The principal coordinates of
distances". `do_pcoa` places the individuals of a `Distances` on components
from the distance of every pair of them, and refuses distances that are not
Euclidean, those that no space has points at; `correct_dists_by_lingoes`
makes them Euclidean by adding one constant to every squared distance.

The comparison with pyNei at commit ef0ca6e, which `pyproject.toml` names,
is on the Kosman distances of `tests/reference/dists/four_alleles.vcf.gz`,
40 individuals, which are Euclidean: the 39 components of popnei against
the first 39 of pyNei, whose 40th is the eigenvalue 0 of the centering. A
component multiplied by -1 is the same component, so pyNei's numbers are
given the sign rule of the spec before they are compared: in each component
the projection of the largest absolute value is made positive.

The worked example is the ten distances of pyNei's `test_pcoa`, of five
individuals, which are not Euclidean. Its numbers are those of R 4.6.1 with
ape 5.8.1 that "How it is verified" of the spec gives and
`tests/reference/pca/small.lingoes.r.*.tsv` holds.

`do_pcoa_from_variants` is the PCoA of the Kosman distances of a `Variants`
in one pass. It is checked on `tests/reference/dists/panel.vcf.gz`, 200
individuals and 1200 variants, whose distances are not Euclidean: with
`correct_by_lingoes` against pyNei's `do_pcoa` of popnei's
`correct_dists_by_lingoes` of popnei's `calc_pairwise_kosman_dists` of the
same file, since pyNei has no correction, and against the literals of R
that the spec gives; without it, refused.
"""

from pathlib import Path

import numpy
import pytest
from popnei import (
    Distances,
    FilteringStats,
    LingoesCorrection,
    PassStats,
    PCoAResult,
    _core,
    calc_pairwise_kosman_dists,
    correct_dists_by_lingoes,
    do_pcoa,
    do_pcoa_from_variants,
    open_vcf,
)
from pynei.dists import Distances as PyneiDistances
from pynei.pca import do_pcoa as pynei_do_pcoa

FOUR_ALLELES_VCF = Path(__file__).parent / "reference" / "dists" / "four_alleles.vcf.gz"
FOUR_ALLELES_NUM_INDIVIDUALS = 40
PANEL_VCF = Path(__file__).parent / "reference" / "dists" / "panel.vcf.gz"

# The tolerance of "How it is verified" of the spec, as an absolute
# difference. The largest number compared is a percentage of the worked
# example, 77.13, where 1e-9 is the eleventh significant digit, and the
# spec's literals carry fifteen. The review of the spec found popnei and
# pyNei 2e-14 apart on the Kosman distances of `four_alleles.vcf.gz`.
TOLERANCE = 1e-9

# The ten distances of pyNei's `test_pcoa`, of the individuals i1 to i5, in
# the order (i1, i2), (i1, i3), ..., (i4, i5).
TEN_DISTS = [0.2, 0.3, 0.9, 0.9, 0.1, 0.8, 0.7, 0.7, 0.8, 0.2]
FIVE_NAMES = ("i1", "i2", "i3", "i4", "i5")

# What Lingoes' correction of the ten gives, from R's ape: the constant c
# added twice to every squared distance, the share of the negative
# eigenvalues of the distances given, and sqrt(0.2² + 2c), the first
# corrected distance.
LINGOES_CONSTANT = 0.0640069399611263
NEGATIVE_EIGENVALUES_PERCENT = 7.88262807403034
FIRST_CORRECTED_DIST = 0.409894962060102

# The principal coordinates of the corrected ten, from R's ape.
CORRECTED_PROJECTIONS = [
    [-0.431869046368213, -0.0415086068851603, 0.224881556185394],
    [-0.283479006142767, -0.158680403878776, -0.138801060767122],
    [-0.269028184151739, 0.212468396790800, -0.131478395434635],
    [0.492920681079785, 0.195146288470815, 0.0612591236372794],
    [0.491455555582935, -0.207425674497679, -0.0158612236209160],
]
CORRECTED_PERCENT = [77.1278402980914, 14.3397713765961, 8.53238832531248]


@pytest.fixture(scope="module")
def kosman_dists() -> Distances:
    """The Kosman distances of the 40 individuals of `four_alleles.vcf.gz`,
    with the names and the counts of their pass."""
    return calc_pairwise_kosman_dists(open_vcf(FOUR_ALLELES_VCF))


def ten_dists(**changed) -> Distances:
    """The ten distances of the worked example, named i1 to i5, with the
    distances at the positions of `changed` replaced by its values."""
    vector = numpy.array(TEN_DISTS)
    for position, value in changed.items():
        vector[int(position.removeprefix("at_"))] = value
    return Distances(vector, names=FIVE_NAMES)


def the_signs_fixed(projections: numpy.ndarray) -> numpy.ndarray:
    """pyNei's projections under the sign rule of the spec: in each
    component the projection of the largest absolute value is made
    positive, the first of them when two have the same absolute value."""
    projections = numpy.array(projections)
    for component in range(projections.shape[1]):
        largest = numpy.argmax(numpy.abs(projections[:, component]))
        if projections[largest, component] < 0:
            projections[:, component] *= -1
    return projections


def test_the_39_components_of_the_kosman_distances_are_pyneis(kosman_dists):
    result = do_pcoa(kosman_dists)
    of_pynei = pynei_do_pcoa(
        PyneiDistances(
            numpy.array(kosman_dists.dist_vector), names=list(kosman_dists.names)
        )
    )
    assert result.projections.shape == (FOUR_ALLELES_NUM_INDIVIDUALS, 39)
    assert of_pynei.projections.shape == (FOUR_ALLELES_NUM_INDIVIDUALS, 40)
    numpy.testing.assert_allclose(
        result.projections.to_numpy(),
        the_signs_fixed(of_pynei.projections.to_numpy()[:, :39]),
        rtol=0,
        atol=TOLERANCE,
    )
    numpy.testing.assert_allclose(
        result.explained_variance_percent.to_numpy(),
        of_pynei.explained_variance_percent.to_numpy()[:39],
        rtol=0,
        atol=TOLERANCE,
    )
    assert list(result.projections.index) == list(of_pynei.projections.index)


def test_the_result_names_the_components_and_the_individuals(kosman_dists):
    result = do_pcoa(kosman_dists)
    names = [f"PC{number:02d}" for number in range(39)]
    assert isinstance(result, PCoAResult)
    assert list(result.projections.columns) == names
    assert list(result.explained_variance_percent.index) == names
    assert list(result.projections.index) == list(kosman_dists.names)
    assert result.explained_variance_percent.sum() == pytest.approx(100, abs=1e-12)
    assert result.explained_variance_percent.iloc[0] == pytest.approx(
        4.33989428824864, abs=TOLERANCE
    )
    assert result.lingoes_constant == 0
    assert result.negative_eigenvalues_percent == 0
    assert result.pass_stats is kosman_dists.pass_stats
    assert result.pass_stats.num_vars > 0


def test_the_ten_distances_are_refused_with_the_correction_named():
    with pytest.raises(ValueError) as refused:
        do_pcoa(ten_dists())
    said = str(refused.value)
    assert "1 of the 5 eigenvalues" in said
    assert "7.88 percent" in said
    assert "`correct_dists_by_lingoes`" in said
    assert "not Euclidean" in said


def test_the_correction_of_the_ten_distances_is_apes():
    given = ten_dists()
    corrected = correct_dists_by_lingoes(given)
    assert isinstance(corrected, LingoesCorrection)
    assert corrected.constant == pytest.approx(LINGOES_CONSTANT, abs=TOLERANCE)
    assert corrected.negative_eigenvalues_percent == pytest.approx(
        NEGATIVE_EIGENVALUES_PERCENT, abs=TOLERANCE
    )
    assert isinstance(corrected.dists, Distances)
    assert corrected.dists.dist_vector[0] == pytest.approx(
        FIRST_CORRECTED_DIST, abs=TOLERANCE
    )
    # Every distance is sqrt(d² + 2c), the first one included.
    numpy.testing.assert_allclose(
        corrected.dists.dist_vector,
        numpy.sqrt(numpy.array(TEN_DISTS) ** 2 + 2 * LINGOES_CONSTANT),
        rtol=0,
        atol=TOLERANCE,
    )
    assert corrected.dists.names == FIVE_NAMES
    assert corrected.dists.pass_stats is None
    # The distances given are as they were.
    assert list(given.dist_vector) == TEN_DISTS


def test_the_corrected_ten_distances_give_the_table_of_ape():
    result = do_pcoa(correct_dists_by_lingoes(ten_dists()).dists)
    assert list(result.projections.index) == list(FIVE_NAMES)
    assert list(result.projections.columns) == ["PC0", "PC1", "PC2"]
    numpy.testing.assert_allclose(
        result.projections.to_numpy(), CORRECTED_PROJECTIONS, rtol=0, atol=TOLERANCE
    )
    numpy.testing.assert_allclose(
        result.explained_variance_percent.to_numpy(),
        CORRECTED_PERCENT,
        rtol=0,
        atol=TOLERANCE,
    )
    # `do_pcoa` corrected nothing itself.
    assert result.lingoes_constant == 0
    assert result.negative_eigenvalues_percent == 0
    assert result.pass_stats is None


def test_the_correction_keeps_the_pass_stats_of_the_distances(kosman_dists):
    given = Distances(
        numpy.array(TEN_DISTS), names=FIVE_NAMES, pass_stats=kosman_dists.pass_stats
    )
    corrected = correct_dists_by_lingoes(given)
    assert corrected.dists.pass_stats is kosman_dists.pass_stats
    assert do_pcoa(corrected.dists).pass_stats is kosman_dists.pass_stats


def test_the_correction_of_euclidean_distances_changes_nothing(kosman_dists):
    corrected = correct_dists_by_lingoes(kosman_dists)
    assert corrected.constant == 0
    assert corrected.negative_eigenvalues_percent == 0
    assert numpy.array_equal(corrected.dists.dist_vector, kosman_dists.dist_vector)
    assert corrected.dists.names == kosman_dists.names
    assert corrected.dists.pass_stats is kosman_dists.pass_stats


def test_distances_built_without_names_are_named_by_their_positions():
    result = do_pcoa(correct_dists_by_lingoes(Distances(TEN_DISTS)).dists)
    assert list(result.projections.index) == [0, 1, 2, 3, 4]


@pytest.mark.parametrize("function", [do_pcoa, correct_dists_by_lingoes])
def test_the_pairs_with_no_distance_are_refused_by_their_names(function):
    # The pairs (i1, i3), at 1, (i3, i4), at 7, and (i3, i5), at 8, have no
    # distance, so i3 is in all three and i1, i4 and i5 in one each.
    dists = ten_dists(at_1=numpy.nan, at_7=numpy.nan, at_8=numpy.nan)
    with pytest.raises(ValueError) as refused:
        function(dists)
    said = str(refused.value)
    assert said.startswith(
        "3 of the 10 pairs of individuals have no distance, the first of them "
        "'i1' and 'i3', and 'i3' is in 3 of them;"
    )
    assert "has to be given a distance" in said
    assert "position" not in said


def test_one_pair_with_no_distance_has_none():
    # The pair (i2, i3), at 4, alone.
    said = refusal_of(ten_dists(at_4=numpy.nan))
    assert said.startswith(
        "1 of the 10 pairs of individuals has no distance, the first of them "
        "'i2' and 'i3', and 'i2' is in 1 of them;"
    ), said


def test_the_individual_named_is_the_first_of_two_in_as_many_pairs():
    # (i2, i4), at 5, and (i4, i5), at 9: i4 is in both, the others in one.
    # (i1, i2), at 0, and (i3, i5), at 8: every individual of the four is
    # in one, and i1 is the first of them.
    assert "'i4' is in 2 of them" in refusal_of(
        ten_dists(at_5=numpy.nan, at_9=numpy.nan)
    )
    assert "'i1' is in 1 of them" in refusal_of(
        ten_dists(at_0=numpy.nan, at_8=numpy.nan)
    )


def refusal_of(dists: Distances) -> str:
    """What `do_pcoa` says when it refuses `dists`."""
    with pytest.raises(ValueError) as refused:
        do_pcoa(dists)
    return str(refused.value)


def test_a_pair_with_no_distance_is_refused_before_the_matrix_is_decomposed():
    # The ten are not Euclidean, and the pair with no distance is what the
    # user hears of.
    assert "has no distance" in refusal_of(ten_dists(at_2=numpy.nan))


@pytest.mark.parametrize("value", [-0.1, numpy.inf])
@pytest.mark.parametrize("function", [do_pcoa, correct_dists_by_lingoes])
def test_a_negative_or_infinite_distance_is_refused_by_the_names_of_its_pair(
    function, value
):
    # The distance at 3 is that of the pair (i1, i5).
    with pytest.raises(ValueError, match="finite and 0 or above") as refused:
        function(ten_dists(at_3=value))
    said = str(refused.value)
    assert said.startswith(f"the distance of 'i1' and 'i5' is {value},")
    assert "F_ST or f_2" in said
    assert "position" not in said


@pytest.mark.parametrize("function", [do_pcoa, correct_dists_by_lingoes])
def test_distances_that_are_all_zero_are_refused(function):
    with pytest.raises(ValueError, match="every distance is 0"):
        function(Distances(numpy.zeros(10), names=FIVE_NAMES))


@pytest.mark.parametrize("function", [do_pcoa, correct_dists_by_lingoes])
def test_one_individual_is_refused(function):
    with pytest.raises(ValueError, match="there is 1 individual, and"):
        function(Distances(numpy.zeros(0), names=("i1",)))


def test_a_constant_beyond_a_float64_is_refused():
    with pytest.raises(ValueError, match="divide the distances"):
        correct_dists_by_lingoes(
            Distances(numpy.array(TEN_DISTS) * 1e200, names=FIVE_NAMES)
        )


@pytest.mark.parametrize("function", [do_pcoa, correct_dists_by_lingoes])
def test_what_is_not_a_distances_is_refused(function):
    square = Distances(TEN_DISTS, names=FIVE_NAMES).square_dists
    with pytest.raises(TypeError, match="Distances.from_square_dists"):
        function(square)


@pytest.mark.parametrize("function", [_core.pcoa, _core.correct_dists_by_lingoes])
def test_a_vector_of_another_length_is_a_defect_of_popnei(function):
    # Only a caller of the private module reaches it: a `Distances` refuses
    # a vector whose length is the pairs of no number of individuals.
    with pytest.raises(RuntimeError, match="holds 10 distances"):
        function(numpy.array(TEN_DISTS), 6)


def test_what_a_user_is_told_to_do_of_a_pair_with_no_distance_is_the_cores():
    # The core writes what to do, which depends on whether the distances were
    # given or came from the variants, and the package puts in the names
    # alone: the exception of the core carries that text as its last
    # argument, and the message a user reads ends with it.
    vector = numpy.array(TEN_DISTS)
    vector[4] = numpy.nan
    with pytest.raises(_core.PcoaPairsWithNoDistance) as of_the_core:
        _core.pcoa(vector, 5)
    remedy = of_the_core.value.args[-1]
    assert remedy.startswith("a principal coordinate analysis places every")
    assert refusal_of(ten_dists(at_4=numpy.nan)).endswith("; " + remedy)


def test_what_a_distance_has_to_be_is_the_cores():
    vector = numpy.array(TEN_DISTS)
    vector[4] = -1.0
    with pytest.raises(_core.PcoaDistanceOutOfRange) as of_the_core:
        _core.pcoa(vector, 5)
    what_it_has_to_be = of_the_core.value.args[-1]
    assert what_it_has_to_be.startswith("a principal coordinate analysis needs")
    assert refusal_of(ten_dists(at_4=-1.0)).endswith(", and " + what_it_has_to_be)


def test_the_corrected_distances_have_no_standard_errors():
    # The standard errors of the distances given are of those distances, and
    # the correction changes every one of them, so the corrected distances
    # carry none.
    dists = Distances(
        numpy.array(TEN_DISTS),
        names=FIVE_NAMES,
        standard_errors=numpy.full(10, 0.01),
    )
    corrected = correct_dists_by_lingoes(dists).dists
    assert corrected.standard_errors is None


# The panel corrected inside the analysis, from R's ape: 198 components, the
# constant c, the share of the negative eigenvalues before the correction,
# the projections of three individuals on the first three components, and
# the first three percentages.
PANEL_NUM_COMPS = 198
PANEL_LINGOES_CONSTANT = 0.014182298472042
PANEL_NEGATIVE_EIGENVALUES_PERCENT = 2.98343616373556
PANEL_PROJECTIONS = {
    "s000": [0.0131009923566961, 0.103593034190948, -0.0461610570616341],
    "s001": [0.0188069490905941, 0.102449570787996, -0.0506243303501242],
    "s199": [-0.0728375733863349, -0.0163081366424616, 0.0108102147666687],
}
PANEL_PERCENT = [9.62407114041929, 6.65101405201371, 1.73580899943624]


@pytest.fixture(scope="module")
def panel_corrected() -> PCoAResult:
    """The PCoA of the variants of the panel with Lingoes' correction."""
    return do_pcoa_from_variants(open_vcf(PANEL_VCF), correct_by_lingoes=True)


def test_the_panel_corrected_is_pyneis_pcoa_of_the_corrected_distances(
    panel_corrected,
):
    # pyNei has no correction, so it is given the distances popnei corrects,
    # and it gives all 200 components, of which the last two are of the two
    # eigenvalues 0: that of the centering and the most negative one before
    # the correction.
    corrected = correct_dists_by_lingoes(
        calc_pairwise_kosman_dists(open_vcf(PANEL_VCF))
    ).dists
    of_pynei = pynei_do_pcoa(
        PyneiDistances(numpy.array(corrected.dist_vector), names=list(corrected.names))
    )
    assert panel_corrected.projections.shape == (200, PANEL_NUM_COMPS)
    assert of_pynei.projections.shape == (200, 200)
    numpy.testing.assert_allclose(
        panel_corrected.projections.to_numpy(),
        the_signs_fixed(of_pynei.projections.to_numpy()[:, :PANEL_NUM_COMPS]),
        rtol=0,
        atol=TOLERANCE,
    )
    numpy.testing.assert_allclose(
        panel_corrected.explained_variance_percent.to_numpy(),
        of_pynei.explained_variance_percent.to_numpy()[:PANEL_NUM_COMPS],
        rtol=0,
        atol=TOLERANCE,
    )
    assert list(panel_corrected.projections.index) == list(of_pynei.projections.index)


def test_the_panel_corrected_gives_the_numbers_of_ape(panel_corrected):
    names = [f"PC{number:03d}" for number in range(PANEL_NUM_COMPS)]
    assert isinstance(panel_corrected, PCoAResult)
    assert list(panel_corrected.projections.columns) == names
    assert list(panel_corrected.explained_variance_percent.index) == names
    assert list(panel_corrected.projections.index) == list(
        open_vcf(PANEL_VCF).individuals
    )
    assert panel_corrected.lingoes_constant == pytest.approx(
        PANEL_LINGOES_CONSTANT, abs=TOLERANCE
    )
    assert panel_corrected.negative_eigenvalues_percent == pytest.approx(
        PANEL_NEGATIVE_EIGENVALUES_PERCENT, abs=TOLERANCE
    )
    for individual, expected in PANEL_PROJECTIONS.items():
        numpy.testing.assert_allclose(
            panel_corrected.projections.loc[individual].to_numpy()[:3],
            expected,
            rtol=0,
            atol=TOLERANCE,
        )
    numpy.testing.assert_allclose(
        panel_corrected.explained_variance_percent.to_numpy()[:3],
        PANEL_PERCENT,
        rtol=0,
        atol=TOLERANCE,
    )
    assert panel_corrected.explained_variance_percent.sum() == pytest.approx(
        100, abs=1e-12
    )


def test_the_panel_counts_the_1200_variants_of_its_pass(panel_corrected):
    assert panel_corrected.pass_stats == PassStats(num_vars=1200, filtering={})


def test_the_panel_counts_the_variants_the_maf_filter_kept():
    # pyNei's `filter_by_maf` at 0.7 keeps 566 of the 1200 variants of the
    # panel as well.
    variants = open_vcf(PANEL_VCF)
    variants.filter_by_maf(0.7)
    result = do_pcoa_from_variants(variants, correct_by_lingoes=True)
    assert result.pass_stats == PassStats(
        num_vars=566,
        filtering={"maf": FilteringStats(vars_processed=1200, vars_kept=566)},
    )


def test_the_panel_is_refused_without_the_correction_named():
    with pytest.raises(ValueError) as refused:
        do_pcoa_from_variants(open_vcf(PANEL_VCF))
    said = str(refused.value)
    assert said.startswith(f"{PANEL_VCF}: 44 of the 200 eigenvalues"), said
    assert "2.98 percent" in said
    assert "not Euclidean" in said
    assert "`correct_by_lingoes`" in said
    assert "correct_dists_by_lingoes" not in said


@pytest.mark.parametrize("correct_by_lingoes", [False, True])
def test_the_pairs_called_together_at_too_few_variants_are_refused_by_their_names(
    correct_by_lingoes,
):
    # With 1105 variants needed, 35 of the 19900 pairs of the panel have no
    # distance, the first of them s001 and s082, and s082 is in 17 of them;
    # they are refused before the matrix is decomposed, corrected or not.
    with pytest.raises(ValueError) as refused:
        do_pcoa_from_variants(
            open_vcf(PANEL_VCF),
            min_num_snps=1105,
            correct_by_lingoes=correct_by_lingoes,
        )
    said = str(refused.value)
    assert said.startswith(
        f"{PANEL_VCF}: 35 of the 19900 pairs of individuals have no distance, "
        f"the first of them 's001' and 's082', and 's082' is in 17 of them; "
    ), said
    assert "`min_num_snps`" in said
    assert "`filter_individuals`" in said
    assert "position" not in said


def test_the_pairs_with_no_distance_are_named_among_the_individuals_kept():
    # Without s000 the positions of the core are one lower than in the
    # source, and the names are still s001, s082 and s082: 17 of the 35 pairs
    # are of s082, and none of them of s000, whose pairs all have a distance.
    variants = open_vcf(PANEL_VCF)
    variants.filter_individuals([f"s{number:03d}" for number in range(1, 200)])
    with pytest.raises(ValueError) as refused:
        do_pcoa_from_variants(variants, min_num_snps=1105)
    said = str(refused.value)
    assert "35 of the 19701 pairs" in said, said
    assert "the first of them 's001' and 's082', and 's082' is in 17" in said


def test_min_num_snps_is_refused_as_the_kosman_distances_refuse_it():
    with pytest.raises(TypeError, match="`min_num_snps` is 1.5"):
        do_pcoa_from_variants(open_vcf(PANEL_VCF), min_num_snps=1.5)
    with pytest.raises(ValueError, match="`min_num_snps` is -1"):
        do_pcoa_from_variants(open_vcf(PANEL_VCF), min_num_snps=-1)


def test_what_is_not_a_variants_is_refused():
    with pytest.raises(TypeError, match="`do_pcoa_from_variants`"):
        do_pcoa_from_variants(PANEL_VCF)


def test_one_individual_is_refused_before_the_pass():
    variants = open_vcf(PANEL_VCF)
    variants.filter_individuals(["s000"])
    with pytest.raises(ValueError, match="there is 1 individual, and") as refused:
        do_pcoa_from_variants(variants)
    assert str(refused.value).startswith(f"{PANEL_VCF}: ")


def test_the_correction_of_the_variants_is_off_by_default_as_the_core_says():
    # The owner decided on 27 September 2026 that distances that are not
    # Euclidean are refused unless the user asks for the correction, and the
    # default is the core's constant.
    import inspect

    default = (
        inspect.signature(do_pcoa_from_variants)
        .parameters["correct_by_lingoes"]
        .default
    )
    assert default is _core.DEFAULT_CORRECT_BY_LINGOES
    assert _core.DEFAULT_CORRECT_BY_LINGOES is False


@pytest.mark.parametrize("given", [1, "yes", None, numpy.int64(1)])
def test_a_correct_by_lingoes_that_is_not_a_bool_is_refused_by_its_name(given):
    with pytest.raises(TypeError, match="`correct_by_lingoes` is True or False"):
        do_pcoa_from_variants(open_vcf(PANEL_VCF), correct_by_lingoes=given)


def test_a_correct_by_lingoes_of_numpy_is_taken():
    result = do_pcoa_from_variants(
        open_vcf(PANEL_VCF), correct_by_lingoes=numpy.bool_(True)
    )
    assert result.projections.shape == (200, PANEL_NUM_COMPS)


# Three individuals of which the first two are never called at one variant,
# so their pair has no distance whatever `min_num_snps` is.
NEVER_CALLED_TOGETHER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ti0\ti1\ti2\n"
    "1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t./.\t1/1\n"
    "1\t2\t.\tA\tC\t.\t.\t.\tGT\t./.\t0/1\t0/0\n"
    "1\t3\t.\tA\tC\t.\t.\t.\tGT\t1/1\t./.\t0/1\n"
    "1\t4\t.\tA\tC\t.\t.\t.\tGT\t./.\t1/1\t0/1\n"
)


@pytest.mark.parametrize("min_num_snps", [None, 0, 1])
def test_a_pair_called_together_at_no_variant_is_not_sent_to_min_num_snps(
    tmp_path, min_num_snps
):
    vcf = tmp_path / "never_called_together.vcf"
    vcf.write_text(NEVER_CALLED_TOGETHER)
    with pytest.raises(ValueError) as refused:
        do_pcoa_from_variants(open_vcf(vcf), min_num_snps=min_num_snps)
    said = str(refused.value)
    assert "'i0' and 'i1'" in said, said
    assert "called together at no variant" in said, said
    assert "min_num_snps" not in said, said
