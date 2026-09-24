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
    _core,
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


def _reference_spectrum() -> dict[str, list[str]]:
    """`panel_folded_sfs_dadi.tsv`, as the column of each population: one
    value per count of the rarer allele, 0 first.

    Its rows are the counts of the rarer allele and not the populations, so
    `_reference` above cannot read it.
    """
    lines = (DIVERSITY_REFERENCE_DIR / "panel_folded_sfs_dadi.tsv").read_text()
    lines = lines.splitlines()
    pops = lines[0].split("\t")[1:]
    of_each_pop: dict[str, list[str]] = {pop: [] for pop in pops}
    for rarer_allele, line in enumerate(lines[1:]):
        values = line.split("\t")
        assert int(values[0]) == rarer_allele
        for pop, value in zip(pops, values[1:], strict=True):
            of_each_pop[pop].append(value)
    return of_each_pop


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


def test_a_call_that_names_no_statistic_gives_the_four_that_need_no_draw() -> None:
    """`calc_pop_diversity(variants)`, with no other argument, gives the four
    statistics that need no draw and no spectrum.

    Those four are what `PopDiversityStat.WITHOUT_A_DRAW` holds. The bins of
    the spectrum are counts of the rarer allele in a draw of
    `num_called_alleles`, so while the default was all five the plainest call
    of the module refused itself and its message named an argument the user
    had not written.
    """
    assert PopDiversityStat.WITHOUT_A_DRAW == (
        PopDiversityStat.NUM_ALLELES,
        PopDiversityStat.PRIVATE_ALLELES,
        PopDiversityStat.VARIABLE_VARS_RATIO,
        PopDiversityStat.FIS,
    )

    diversity = calc_pop_diversity(_panel())

    assert diversity.pops == ("pop",)
    assert diversity.num_alleles.loc["pop", "total"] > 0
    assert diversity.private_alleles.loc["pop", "total"] > 0
    assert diversity.variable_vars_ratio.loc["pop", "total"] > 0
    assert not math.isnan(diversity.fis["pop"])
    # The one statistic that needs a draw is not in the default, and naming
    # it without a draw is refused as it was.
    assert diversity.folded_sfs is None
    assert diversity.num_vars.loc["pop", "with_data"] == PANEL_NUM_VARS


def test_the_default_stats_are_the_ones_of_the_core() -> None:
    """The four statistics that need no draw are named in the Rust core, and
    the package reads them from it, so that one added there is in the default
    of both packages with nothing written in either."""
    assert PopDiversityStat.WITHOUT_A_DRAW == tuple(
        PopDiversityStat(name) for name in _core.diversity_stats_without_a_draw()
    )
    assert PopDiversityStat.FOLDED_SFS not in PopDiversityStat.WITHOUT_A_DRAW


def test_the_statistics_are_given_by_name_and_not_after_the_populations() -> None:
    """`stats` and the two arguments after it are keyword only, so that the
    argument after the variants cannot mean the populations in one function of
    popnei and the statistics in another. `variants` and `pops` stay
    positional."""
    with pytest.raises(TypeError, match="positional"):
        calc_pop_diversity(_panel(), PANEL_POPS, WITH_NO_SPECTRUM)

    diversity = calc_pop_diversity(_panel(), PANEL_POPS, stats=WITH_NO_SPECTRUM)

    assert diversity.pops == PANEL_POP_NAMES


def test_the_draw_and_the_spectrum_are_not_built_yet() -> None:
    """A `num_called_alleles` and a `stats` that names the spectrum are
    refused while the pass computes neither, and the message names the work
    package that builds them.

    Without the refusal both give an answer that is wrong and reads as an
    answer: NaN in the three `in_draw` columns, which this module keeps for a
    draw above every population's called alleles, and `None` in `folded_sfs`,
    which it keeps for a statistic nobody asked for.
    """
    with pytest.raises(NotImplementedError, match="work package 3"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=20)
    with pytest.raises(NotImplementedError, match="work package 3"):
        _of_the_panel(stats=(PopDiversityStat.FOLDED_SFS,))
    with pytest.raises(NotImplementedError, match="folded site frequency spectrum"):
        _of_the_panel(stats=(PopDiversityStat.FIS, PopDiversityStat.FOLDED_SFS))


def test_a_num_called_alleles_that_is_no_draw_is_refused_as_a_wrong_argument() -> None:
    """A value that is no draw at all is wrong whatever popnei computes, so it
    is refused as one and not as a statistic that is not built: a whole number
    is what the argument takes, and 0 and 1 are draws no standardized value
    can be taken over."""
    with pytest.raises(TypeError, match="num_called_alleles"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=2.5)
    with pytest.raises(TypeError, match="num_called_alleles"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=True)
    with pytest.raises(ValueError, match="num_called_alleles"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=-1)
    with pytest.raises(ValueError, match="num_called_alleles"):
        _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=1)


@pytest.mark.xfail(
    strict=True,
    reason="the draw is task 3.5 of `docs/plans/diversity.md`, and until it is "
    "built `calc_pop_diversity` refuses a `num_called_alleles`. The day the pass "
    "fills the three `in_draw` columns and the two counts of the variants in a "
    "draw this passes, `strict` turns that into a failure of the suite, and "
    "whoever built it has to take this marker off and leave the assertions.",
)
def test_the_standardized_values_of_the_panel_are_vegans() -> None:
    """The alleles a draw of 20 is expected to show and the chance that such a
    draw varies, averaged over the variants of each population of the panel,
    against `vegan`'s numbers stored in `tests/reference/diversity/`; and the
    totals, which no draw changes.

    The standardized private alleles are not here: no program outside popnei
    computes one, and "How it is verified" of that item says what checks
    them.
    """
    with_no_draw = _of_the_panel(stats=WITH_NO_SPECTRUM)
    of_a_draw_of_20 = _of_the_panel(stats=WITH_NO_SPECTRUM, num_called_alleles=20)
    of_the_alleles = _reference("panel_num_alleles.tsv")
    of_the_variable = _reference("panel_variable_vars.tsv")

    assert of_a_draw_of_20.num_vars_every_pop_in_draw == PANEL_NUM_VARS
    for pop in PANEL_POP_NAMES:
        assert of_a_draw_of_20.num_vars.loc[pop, "in_draw"] == int(
            of_the_alleles[pop]["num_vars_in_draw"]
        )
        assert _the_same_number(
            of_a_draw_of_20.num_alleles.loc[pop, "in_draw"],
            of_the_alleles[pop]["in_draw_vegan"],
        ), f"the alleles a draw of 20 shows in {pop}"
        assert _the_same_number(
            of_a_draw_of_20.variable_vars_ratio.loc[pop, "in_draw"],
            of_the_variable[pop]["in_draw_vegan_minus_one"],
        ), f"the chance that a draw of 20 varies in {pop}"
        # A draw changes the standardized values alone.
        for statistic, column in (
            ("num_alleles", "total"),
            ("private_alleles", "total"),
            ("variable_vars_ratio", "total"),
        ):
            assert (
                getattr(of_a_draw_of_20, statistic).loc[pop, column]
                == getattr(with_no_draw, statistic).loc[pop, column]
            ), f"the {column} of the {statistic} of {pop}"
        assert of_a_draw_of_20.fis[pop] == with_no_draw.fis[pop]


@pytest.mark.xfail(
    strict=True,
    reason="the folded site frequency spectrum is task 3.5 of "
    "`docs/plans/diversity.md`, and until it is built `calc_pop_diversity` "
    "refuses a `stats` that names it. The day the pass fills it this passes, "
    "`strict` turns that into a failure of the suite, and whoever built it has "
    "to take this marker off and leave the assertions.",
)
def test_the_folded_spectrum_of_the_panel_is_dadis() -> None:
    """The variants of each population of the panel expected to show each
    count of their rarer allele in a draw of 20, against `dadi`'s numbers
    stored in `tests/reference/diversity/panel_folded_sfs_dadi.tsv`: eleven
    rows, the counts 0 to 10, each column summing to the 1200 variants that
    counted."""
    diversity = _of_the_panel(
        stats=(PopDiversityStat.FOLDED_SFS,), num_called_alleles=20
    )
    of_dadi = _reference_spectrum()

    assert list(diversity.folded_sfs.index) == list(range(11))
    for pop in PANEL_POP_NAMES:
        for rarer_allele, theirs in enumerate(of_dadi[pop]):
            assert _the_same_number(
                diversity.folded_sfs.loc[rarer_allele, pop], theirs
            ), f"the variants of {pop} with {rarer_allele} copies of the rarer allele"
        assert _the_same_number(diversity.folded_sfs[pop].sum(), PANEL_NUM_VARS), (
            f"the spectrum of {pop} over its variants"
        )


@pytest.mark.xfail(
    strict=True,
    reason="the spectrum is task 3.5 of `docs/plans/diversity.md`, and until it "
    "is built `calc_pop_diversity` refuses a `stats` that names it with a "
    "`NotImplementedError`, before the core is reached. The day the pass fills "
    "the spectrum this `ValueError` of the core is what a user who gives no "
    "`num_called_alleles` reads again, `strict` turns the pass into a failure of "
    "the suite, and whoever built it has to take this marker off.",
)
def test_the_spectrum_asked_for_with_no_draw_is_refused() -> None:
    """The bins of a folded spectrum are the counts of the rarer allele in a
    draw, so the spectrum needs `num_called_alleles`. A user who names no
    statistic asks for the four that need no draw and does not reach this."""
    with pytest.raises(ValueError, match="num_called_alleles"):
        _of_the_panel(stats=(PopDiversityStat.FOLDED_SFS,))
    with pytest.raises(ValueError, match="folded site frequency spectrum"):
        _of_the_panel(stats=(PopDiversityStat.FIS, PopDiversityStat.FOLDED_SFS))


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
    # Every population's mean private alleles is over the variants every
    # population counted, which is 0 here, so all three are NaN. `p1` and `p2`
    # counted variants of their own, so they are what tells that divisor from
    # the variants of the population.
    for pop in PANEL_POP_NAMES:
        assert math.isnan(diversity.private_alleles.loc[pop, "mean"]), (
            f"the mean private alleles of {pop}"
        )
    for pop in ("p1", "p2"):
        assert diversity.num_vars.loc[pop, "with_data"] > 0
        assert diversity.num_alleles.loc[pop, "total"] > 0
        assert not math.isnan(diversity.fis[pop])


def test_the_mean_private_alleles_are_over_the_variants_every_population_counted(
    write_vcf,
) -> None:
    """Two variants of three individuals in which one population misses the
    second, so that the variants every population counted are fewer than the
    variants of `pop1` and the two divisors give different numbers.

    `pop1` is `ind1` and `pop2` is `ind2`; `ind3` is in no population and
    takes no part. At the first variant `ind1` is `0/1` and `ind2` is `2/2`,
    so `pop1` called the alleles 0 and 1 and `pop2` called 2, and each of the
    three is private to the population that called it. At the second `ind1` is
    `0/0` and `ind2` has no genotype, so the variant counts for `pop1` alone
    and is out of the private alleles of both.

    `pop1` therefore has 2 variants with data, 2 private alleles and 1 variant
    every population counted, which makes its mean private alleles 2. Over its
    own 2 variants it would be 1, and over the same divisor its 3 alleles
    called and its 1 variable variant would be 3 and 1 rather than 1.5 and
    0.5.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT,G\t.\tPASS\t.\tGT\t0/1\t2/2\t0/0",
            "chr1\t20\t.\tA\tT,G\t.\tPASS\t.\tGT\t0/0\t./.\t0/0",
        ]
    )

    diversity = calc_pop_diversity(
        open_vcf(path),
        pops={"pop1": ["ind1"], "pop2": ["ind2"]},
        stats=WITH_NO_SPECTRUM,
        min_num_individuals=1,
    )

    assert diversity.num_vars.loc["pop1", "with_data"] == 2
    assert diversity.num_vars.loc["pop2", "with_data"] == 1
    assert diversity.num_vars_every_pop == 1
    assert diversity.private_alleles.loc["pop1", "total"] == 2
    assert diversity.private_alleles.loc["pop1", "mean"] == 2.0
    assert diversity.private_alleles.loc["pop2", "total"] == 1
    assert diversity.private_alleles.loc["pop2", "mean"] == 1.0
    assert diversity.num_alleles.loc["pop1", "total"] == 3
    assert diversity.num_alleles.loc["pop1", "mean"] == 1.5
    assert diversity.variable_vars_ratio.loc["pop1", "total"] == 1
    assert diversity.variable_vars_ratio.loc["pop1", "ratio"] == 0.5


def test_a_pass_that_gives_no_variant_is_refused(write_vcf) -> None:
    """Every count of a population is over the variants of the pass, so a pass
    with none is an error, and the message says whether the source held none or
    the steps kept none, with what each filter counted."""
    variants = _panel()
    # No variant with a called allele has a major allele frequency of 0, and
    # one without a called allele is not kept either.
    variants.filter_by_maf(0)

    with pytest.raises(ValueError, match="the pass gave no variant") as refusal:
        calc_pop_diversity(variants, pops=PANEL_POPS, stats=WITH_NO_SPECTRUM)
    assert "the `maf` filter was given 1200 and kept 0" in str(refusal.value)

    with pytest.raises(ValueError, match="its source holds none"):
        calc_pop_diversity(open_vcf(write_vcf([])), stats=WITH_NO_SPECTRUM)


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
    for a statistic nobody asked for, and asking for none is refused by the
    Rust core, which lists the five names a user can write."""
    with pytest.raises(TypeError, match="PopDiversityStat"):
        _of_the_panel(stats=("num_alleles",))
    with pytest.raises(TypeError, match="PopDiversityStat"):
        _of_the_panel(stats="fis")
    with pytest.raises(ValueError, match="`stats` names no statistic") as refusal:
        _of_the_panel(stats=())
    assert "reads every variant of the source for nothing" in str(refusal.value)


def test_a_name_that_is_of_no_statistic_is_refused_and_names_no_file() -> None:
    """What a user who calls `popnei._core` themselves reads.

    The package takes the members of `PopDiversityStat` and nothing else, so a
    name reaches the Rust core only from a caller that went round the package.
    The core lists the five names a user can write, and the message names no
    file: what a user wrote is wrong whatever variants are read.
    """
    variants = _panel()

    with pytest.raises(ValueError, match="`num_allelez` is not one of the") as refusal:
        _core.calc_pop_diversity(
            variants._source, variants._steps, None, ["num_allelez"], None, 20
        )

    message = str(refusal.value)
    assert "num_alleles, private_alleles, variable_vars_ratio, folded_sfs, fis" in (
        message
    )
    assert PANEL.name not in message
    assert getattr(refusal.value, "filename", None) is None


def test_one_statistic_written_on_its_own_is_that_one_statistic() -> None:
    """A member of a `StrEnum` is a string, so one written without its comma
    would be a sequence of its letters."""
    diversity = _of_the_panel(stats=PopDiversityStat.NUM_ALLELES)

    assert diversity.num_alleles is not None
    assert diversity.fis is None
