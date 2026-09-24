"""How much variety each population holds, from Python: the alleles each
population called, the private ones among them, the variants that vary in it
and F_IS, over the panel and over the cases a user can reach.

`docs/specs/diversity.md` has the five statistics and, under "How it is
verified" of each, the program its numbers come from. pyNei has none of the
five, so no test here runs the two libraries on the same input.

The dataset is the panel of `docs/specs/stats.md`,
`tests/reference/stats/panel.vcf.gz`: 1200 biallelic diploid variants of 200
individuals, 3 in 100 genotypes missing whole, in the three populations `p0`,
`p1` and `p2` of 48, 68 and 84 individuals that `panel_pops_bcftools.txt`
beside it holds. At `min_num_individuals` 20 all 1200 variants count for all
three.

The counts the tests assert are read from the files of
`tests/reference/diversity/`, whose `README.md` says which file holds what
and which program wrote it: `adegenet` 2.1.11 counted the alleles and the
variable variants, `poppr` 2.9.8 the private ones, and `scikit-allel` 1.3.13
gave the plain form of F_IS. The unbiased form, which popnei returns, is in
the spec and not in a file: no program outside popnei computes it, and the
values there come from `docs/reports/diversity-method/panel.py`, which works
the five quantities out in Python as the spec defines them.

The draw of a common number of called alleles is not asserted here: the
standardized values and the folded spectrum are work package 3 of
`docs/plans/diversity.md`.
"""

import math
from pathlib import Path

import pytest
from popnei import (
    PerVarStat,
    PopDiversity,
    PopDiversityStat,
    calc_per_var_distribs,
    calc_pop_diversity,
    open_vcf,
)
from popnei.variant import Variants

STATS_REFERENCE_DIR = Path(__file__).parent / "reference" / "stats"
DIVERSITY_REFERENCE_DIR = Path(__file__).parent / "reference" / "diversity"
PANEL = STATS_REFERENCE_DIR / "panel.vcf.gz"

# The variants of the panel, and how many called genotypes a population needs
# at one of them for the variant to count for it, which "How it is verified"
# of every item of the spec ran the reference programs at.
PANEL_NUM_VARS = 1200
PANEL_MIN_NUM_INDIVIDUALS = 20

# What a value of popnei may differ from the number it is compared with by,
# which is what every item of the spec asks for its floats: 1e-12 of the
# value. Every number compared here carries all its digits, the stored ones
# because the files of `tests/reference/diversity/` hold 17 significant
# figures and the three F_IS because the spec gives them at the precision
# they were computed to, and the two sides add the same per variant values in
# different orders, so the last bits differ.
OF_A_REFERENCE_VALUE = 1e-12

# The unbiased F_IS of the three populations of the panel, from "How it is
# verified" of "The inbreeding coefficient F_IS" of the spec. It is the form
# popnei returns: one minus the mean observed heterozygosity over the mean
# unbiased expected one.
#
# The item prints those three to ten decimals for the reader and gives them
# again at the precision they were computed to, which is what is written
# here: a value of ten decimals stands for anything within 5e-11 of itself,
# and 1e-12 of 0.0128 is 1.3e-14, so the shorter form could not be compared
# within the tolerance this file uses everywhere. They come from
# `docs/reports/diversity-method/panel.py`, which computes the five
# quantities of the spec in Python and is where every unbiased number of the
# spec comes from, so this assertion says that popnei's Rust agrees with that
# Python and nothing more, as the standardized private alleles do. What
# checks F_IS against a program outside popnei is the plain form below,
# against `scikit-allel`.
PANEL_FIS = {
    "p0": -0.012758486763377208,
    "p1": -0.018110713076467277,
    "p2": -0.018458583231322434,
}


def _pops_of(path: Path) -> dict[str, list[str]]:
    """The populations of a file of individual name and population name, one
    pair per line, in the order the file has them."""
    pops: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        name, pop = line.split("\t")
        pops.setdefault(pop, []).append(name)
    return pops


PANEL_POPS = _pops_of(STATS_REFERENCE_DIR / "panel_pops_bcftools.txt")

# The three populations in the order the keys of `PANEL_POPS` iterate in,
# which is the order the file names them in and the order every result of a
# call with those `pops` has to be in.
PANEL_POP_NAMES = tuple(PANEL_POPS)


def _reference(name: str) -> dict[str, dict[str, str]]:
    """One file of `tests/reference/diversity/`, as the row of each
    population under the name of its column."""
    lines = (DIVERSITY_REFERENCE_DIR / name).read_text().splitlines()
    columns = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        row = dict(zip(columns, line.split("\t"), strict=True))
        rows[row["pop"]] = row
    return rows


def _panel() -> Variants:
    """The 1200 variants of the panel."""
    return open_vcf(PANEL)


def _of_the_panel(**kwargs) -> PopDiversity:
    """The diversity of the three populations of the panel, with the
    threshold the reference programs were run at unless another is given."""
    kwargs.setdefault("pops", PANEL_POPS)
    kwargs.setdefault("min_num_individuals", PANEL_MIN_NUM_INDIVIDUALS)
    return calc_pop_diversity(_panel(), **kwargs)


def _the_same_number(ours, theirs) -> bool:
    """Whether two values agree within `OF_A_REFERENCE_VALUE` of the value,
    with a missing value on both sides or on neither."""
    ours, theirs = float(ours), float(theirs)
    if math.isnan(ours) or math.isnan(theirs):
        return math.isnan(ours) and math.isnan(theirs)
    return math.isclose(ours, theirs, rel_tol=OF_A_REFERENCE_VALUE, abs_tol=0)


# The four statistics of this work package, which are every one but the
# folded spectrum: that one needs a draw, and the draw is work package 3.
WITH_NO_SPECTRUM = (
    PopDiversityStat.NUM_ALLELES,
    PopDiversityStat.PRIVATE_ALLELES,
    PopDiversityStat.VARIABLE_VARS_RATIO,
    PopDiversityStat.FIS,
)


def test_the_panel_holds_the_alleles_and_the_fis_the_reference_programs_give() -> None:
    """The alleles the three populations of the panel called, the private
    ones among them and the variants that vary in each, against the counts of
    `adegenet` and `poppr` stored in `tests/reference/diversity/`, and F_IS
    against the three values of the spec, which the Python of
    `docs/reports/diversity-method/panel.py` computed."""
    diversity = _of_the_panel(stats=WITH_NO_SPECTRUM)
    of_the_alleles = _reference("panel_num_alleles.tsv")
    of_the_private = _reference("panel_private_alleles.tsv")
    of_the_variable = _reference("panel_variable_vars.tsv")

    assert diversity.pops == PANEL_POP_NAMES
    for pop in PANEL_POP_NAMES:
        fis = PANEL_FIS[pop]
        # The counts are counts and are compared exactly.
        assert diversity.num_vars.loc[pop, "with_data"] == int(
            of_the_alleles[pop]["num_vars_with_data"]
        )
        assert diversity.num_alleles.loc[pop, "total"] == int(
            of_the_alleles[pop]["total_adegenet"]
        )
        assert diversity.private_alleles.loc[pop, "total"] == int(
            of_the_private[pop]["total_poppr"]
        )
        assert diversity.variable_vars_ratio.loc[pop, "total"] == int(
            of_the_variable[pop]["total_adegenet"]
        )
        assert _the_same_number(diversity.fis[pop], fis), (
            f"the F_IS of {pop} is {diversity.fis[pop]!r} and the spec gives {fis}"
        )
    assert diversity.num_vars_every_pop == int(
        of_the_private["p0"]["num_vars_every_pop"]
    )


def test_the_means_and_the_ratios_of_the_panel_are_the_totals_over_its_variants() -> (
    None
):
    """The mean alleles of each population of the panel, its mean private
    alleles and its ratio of variable variants, against the three numbers the
    reference files hold beside each total."""
    diversity = _of_the_panel(stats=WITH_NO_SPECTRUM)
    of_the_alleles = _reference("panel_num_alleles.tsv")
    of_the_private = _reference("panel_private_alleles.tsv")
    of_the_variable = _reference("panel_variable_vars.tsv")

    for pop in PANEL_POP_NAMES:
        assert _the_same_number(
            diversity.num_alleles.loc[pop, "mean"], of_the_alleles[pop]["mean"]
        ), f"the mean alleles of {pop}"
        assert _the_same_number(
            diversity.private_alleles.loc[pop, "mean"], of_the_private[pop]["mean"]
        ), f"the mean private alleles of {pop}"
        assert _the_same_number(
            diversity.variable_vars_ratio.loc[pop, "ratio"],
            of_the_variable[pop]["ratio"],
        ), f"the ratio of variable variants of {pop}"


def test_the_plain_fis_of_the_two_mean_heterozygosities_is_scikit_allels() -> None:
    """One minus the mean observed heterozygosity over the mean plain
    expected one, built from what `calc_per_var_distribs` gives over the same
    populations and the same threshold, against `scikit-allel`.

    No function of popnei gives that form: `calc_pop_diversity` divides by the
    unbiased expected heterozygosity, which corrects each variant by
    c / (c - 1) with c the called alleles there, and `scikit-allel` does not
    correct it. So this checks the half of F_IS that is the two means and the
    ratio of them, and the other half, the correction, is checked against
    pyNei and plink2 in the tests of `docs/specs/stats.md`.
    """
    distribs = calc_per_var_distribs(
        _panel(),
        stats=(PerVarStat.OBS_HET, PerVarStat.EXP_HET),
        pops=PANEL_POPS,
        min_num_individuals=PANEL_MIN_NUM_INDIVIDUALS,
    )
    of_scikit_allel = _reference("panel_fis_plain_allel.tsv")

    for pop in PANEL_POP_NAMES:
        plain = 1 - distribs.obs_het.mean[pop] / distribs.exp_het.mean[pop]
        assert _the_same_number(plain, of_scikit_allel[pop]["fis_plain"]), (
            f"the plain F_IS of {pop}"
        )


def test_one_population_holds_every_allele_it_called_as_a_private_one() -> None:
    """With no `pops` there is one population of every individual, and with
    no other population to hold them every allele it called is private."""
    diversity = calc_pop_diversity(
        _panel(),
        stats=(PopDiversityStat.NUM_ALLELES, PopDiversityStat.PRIVATE_ALLELES),
        min_num_individuals=PANEL_MIN_NUM_INDIVIDUALS,
    )

    assert diversity.pops == ("pop",)
    total = diversity.num_alleles.loc["pop", "total"]
    assert total > 0
    assert diversity.private_alleles.loc["pop", "total"] == total
    assert diversity.private_alleles.loc["pop", "mean"] == pytest.approx(
        diversity.num_alleles.loc["pop", "mean"]
    )
    assert diversity.num_vars_every_pop == diversity.num_vars.loc["pop", "with_data"]


@pytest.mark.parametrize("poly_threshold", [0.0, 0.5, 0.95, 1.0])
def test_the_variable_variants_are_the_ones_the_stats_module_counts(
    poly_threshold: float,
) -> None:
    """The variable variants of each population are the `num_variable` of the
    polymorphism ratio of `calc_per_var_distribs` over the same populations
    and the same threshold, whatever `poly_threshold` is: that count is the
    variants whose major allele frequency is below 1 and does not read the
    threshold."""
    diversity = _of_the_panel(stats=(PopDiversityStat.VARIABLE_VARS_RATIO,))
    distribs = calc_per_var_distribs(
        _panel(),
        stats=(PerVarStat.POLY_VARS_RATIO,),
        pops=PANEL_POPS,
        min_num_individuals=PANEL_MIN_NUM_INDIVIDUALS,
        poly_threshold=poly_threshold,
    )

    for pop in PANEL_POP_NAMES:
        assert (
            diversity.variable_vars_ratio.loc[pop, "total"]
            == distribs.poly_vars_ratio.num_variable[pop]
        ), f"the variable variants of {pop}"


def test_a_statistic_that_was_not_asked_for_has_no_value() -> None:
    """A result holds the statistics that were asked for, and the others are
    `None`; the counts of the variants are there whatever was asked for."""
    diversity = _of_the_panel(stats=(PopDiversityStat.FIS,))

    assert diversity.fis is not None
    assert diversity.num_alleles is None
    assert diversity.private_alleles is None
    assert diversity.variable_vars_ratio is None
    assert diversity.folded_sfs is None
    assert diversity.num_vars.loc["p0", "with_data"] == PANEL_NUM_VARS


def test_the_spectrum_asked_for_with_no_draw_is_refused() -> None:
    """The bins of a folded spectrum are the counts of the rarer allele in a
    draw, so the spectrum needs `num_called_alleles`. It is the default
    `stats`, so a call that asks for everything and gives no draw is refused
    too."""
    with pytest.raises(ValueError, match="num_called_alleles"):
        _of_the_panel(stats=(PopDiversityStat.FOLDED_SFS,))
    with pytest.raises(ValueError, match="folded site frequency spectrum"):
        _of_the_panel()


def test_a_draw_of_fewer_than_two_alleles_is_refused() -> None:
    """A draw of one allele shows one allele whatever the population holds,
    so every standardized value of it would say nothing."""
    with pytest.raises(ValueError, match="num_called_alleles"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=1)


def test_a_population_that_names_no_individual_of_the_pass_is_refused() -> None:
    """`pops` is read as `calc_per_var_distribs` reads it: a name that is not
    an individual of the pass names the population and the name."""
    with pytest.raises(ValueError, match="s999"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, pops={"p0": ["s000", "s999"]})


def test_a_population_short_of_the_threshold_everywhere_counts_no_variant() -> None:
    """At `min_num_individuals` 60 the 48 individuals of `p0` can never reach
    the threshold, so no variant counts for it: its totals are 0, its mean and
    its ratio are NaN and its F_IS is NaN, and the other two populations are
    given as they are. The divisor of the private alleles, the variants every
    population counted, is 0 as well."""
    diversity = _of_the_panel(stats=WITH_NO_SPECTRUM, min_num_individuals=60)

    assert diversity.num_vars.loc["p0", "with_data"] == 0
    assert diversity.num_alleles.loc["p0", "total"] == 0
    assert diversity.variable_vars_ratio.loc["p0", "total"] == 0
    assert math.isnan(diversity.num_alleles.loc["p0", "mean"])
    assert math.isnan(diversity.variable_vars_ratio.loc["p0", "ratio"])
    assert math.isnan(diversity.fis["p0"])
    assert diversity.num_vars_every_pop == 0
    assert math.isnan(diversity.private_alleles.loc["p0", "mean"])
    for pop in ("p1", "p2"):
        assert diversity.num_vars.loc[pop, "with_data"] > 0
        assert diversity.num_alleles.loc[pop, "total"] > 0
        assert not math.isnan(diversity.fis[pop])


def test_the_totals_and_the_fis_do_not_read_the_draw() -> None:
    """A draw changes the standardized values alone: the alleles called, the
    private ones, the variable variants and F_IS are the same numbers with
    `num_called_alleles` 20 as with none."""
    with_no_draw = _of_the_panel(stats=WITH_NO_SPECTRUM)
    of_a_draw_of_20 = _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=20)

    for pop in PANEL_POP_NAMES:
        assert (
            of_a_draw_of_20.num_alleles.loc[pop, "total"]
            == with_no_draw.num_alleles.loc[pop, "total"]
        )
        assert (
            of_a_draw_of_20.private_alleles.loc[pop, "total"]
            == with_no_draw.private_alleles.loc[pop, "total"]
        )
        assert (
            of_a_draw_of_20.variable_vars_ratio.loc[pop, "total"]
            == with_no_draw.variable_vars_ratio.loc[pop, "total"]
        )
        assert of_a_draw_of_20.fis[pop] == with_no_draw.fis[pop]


def test_the_pass_stats_count_the_variants_the_pass_gave() -> None:
    """Every consumer of popnei gives the counts of its pass, and a pass with
    no filter took the 1200 variants of the panel."""
    diversity = _of_the_panel(stats=WITH_NO_SPECTRUM)

    assert diversity.pass_stats.num_vars == PANEL_NUM_VARS
    assert diversity.pass_stats.filtering == {}


def test_variants_that_is_not_a_variants_is_refused() -> None:
    """The path of the VCF is the mistake that is easiest to make, and what
    it gave was an `AttributeError` about an object with no source inside."""
    with pytest.raises(TypeError, match="open_vcf"):
        calc_pop_diversity(PANEL, pops=PANEL_POPS)


def test_stats_takes_the_members_of_the_enumeration_alone() -> None:
    """A name written as a string is refused, so a typo in one cannot pass
    for a statistic nobody asked for, and asking for none is refused."""
    with pytest.raises(TypeError, match="PopDiversityStat"):
        _of_the_panel(stats=("num_alleles",))
    with pytest.raises(TypeError, match="PopDiversityStat"):
        _of_the_panel(stats="fis")
    with pytest.raises(ValueError, match="names no statistic"):
        _of_the_panel(stats=())


def test_one_statistic_written_on_its_own_is_that_one_statistic() -> None:
    """A member of a `StrEnum` is a string, so one written without its comma
    would be a sequence of its letters."""
    diversity = _of_the_panel(stats=PopDiversityStat.NUM_ALLELES)

    assert diversity.num_alleles is not None
    assert diversity.fis is None
