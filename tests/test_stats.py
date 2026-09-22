"""The per variant distributions from Python: the five statistics of one
pass, per population, and what the call refuses.

`docs/specs/stats.md` has them under "The per variant distributions", and
the comparison it asks for is with pyNei at commit ef0ca6e, which
`pyproject.toml` names, on two datasets. The panel is
`tests/reference/stats/panel.vcf.gz`, 1200 biallelic diploid variants of
200 individuals named `s000` to `s199`, 3 in 100 genotypes missing whole,
in the three populations `p0`, `p1` and `p2` of 48, 68 and 84 individuals
that `panel_pops_bcftools.txt` beside it holds. The other is
`tests/reference/vcf/many.vcf`, 500 variants of 50 diploid individuals,
`ind00` to `ind49`, one in ten with three alleles and 257 half called
genotypes among its 25000, in the two populations `popA`, the first 20
individuals, and `popB`, the other 30, that `many_pops.txt` holds; popnei
reads it with `only_passed=False`, because pyNei gives every variant
whatever its FILTER says.

pyNei is run twice on each dataset, once with `unbiased_exp_het=False`
for popnei's `exp_het` and once with it true for popnei's
`unbiased_exp_het`, which are two statistics here and one with a switch
there. Every value is read by the name of its population and not by its
place, because pyNei sorts the populations of the expected heterozygosity
and keeps the order of the keys for the others, so its columns do not line
up with each other.
"""

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy
import pytest
from popnei import (
    PerVarStat,
    PolyVarsStats,
    StatsDistrib,
    calc_per_var_distribs,
    open_vcf,
)
from popnei.variant import Variants
from pynei import calc_per_var_distribs as pynei_calc_per_var_distribs
from pynei import vars_from_vcf
from pynei.diversity import _calc_unbiased_exp_het_per_var
from pynei.utils_pop import _calc_pops_idxs

STATS_REFERENCE_DIR = Path(__file__).parent / "reference" / "stats"
PANEL = STATS_REFERENCE_DIR / "panel.vcf.gz"
MANY = Path(__file__).parent / "reference" / "vcf" / "many.vcf"

# The four statistics that have a distribution, under the name of the field
# of the result that holds each one.
DISTRIBS = ("obs_het", "maf", "exp_het", "unbiased_exp_het")

# The variants and the individuals of the two datasets.
PANEL_NUM_VARS = 1200
MANY_NUM_VARS = 500

# How many called genotypes a population needs at a variant in each dataset
# for the variant to have a value there, which "How it is verified" of the
# pass gives: the default on the panel, whose populations have 40 called
# genotypes or more at every variant, and 5 on `many.vcf`, whose smallest
# population has 15 at every variant.
PANEL_MIN_NUM_INDIVIDUALS = 20
MANY_MIN_NUM_INDIVIDUALS = 5


def _pops_of(path: Path) -> dict[str, list[str]]:
    """The populations of a file of individual name and population name,
    one pair per line, in the order the file has them."""
    pops: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        name, pop = line.split("\t")
        pops.setdefault(pop, []).append(name)
    return pops


PANEL_POPS = _pops_of(STATS_REFERENCE_DIR / "panel_pops_bcftools.txt")
MANY_POPS = _pops_of(STATS_REFERENCE_DIR / "many_pops.txt")


def _panel() -> Variants:
    """The 1200 variants of the panel."""
    return open_vcf(PANEL)


def _many() -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which is what pyNei reads."""
    return open_vcf(MANY, only_passed=False)


def _the_same_number(ours, theirs) -> bool:
    """Whether two means or two ratios agree within 1e-12 relative, with a
    missing value on both sides or on neither.

    numpy and the Rust loop add the variants of a population in different
    orders, so the last bits of a mean differ.
    """
    ours, theirs = float(ours), float(theirs)
    if math.isnan(ours) or math.isnan(theirs):
        return math.isnan(ours) and math.isnan(theirs)
    # The relative tolerance alone, which "How it is verified" of the pass
    # asks for: an absolute one of 1e-12 is the wider of the two for a value
    # near 0, and it would take any two means below it for equal.
    return math.isclose(ours, theirs, rel_tol=1e-12, abs_tol=0)


def _on_each_bin_edge(path: Path, pops, min_num_individuals: int, edges):
    """How many variants of each population have an unbiased expected
    heterozygosity that lies on each edge of the histogram, from pyNei's own
    per variant values: one count per edge, per population.

    The two libraries reach that value by different arithmetic, popnei by
    multiplying the factors of each of its terms one over another and pyNei
    by multiplying the plain value by `c / (c - 1)`, so the two differ in
    the last bit and such a variant is counted in either of the two bins
    that share the edge, which "How it is verified" of the pass gives with
    the two variants of these datasets that it happens to.
    """
    theirs = vars_from_vcf(path)
    of_each_pop = _calc_pops_idxs(pops, theirs.samples)
    edges = numpy.asarray(edges)
    on_each_edge = {pop: numpy.zeros(len(edges), dtype=int) for pop in of_each_pop}
    for chunk in theirs.iter_vars_chunks():
        values = _calc_unbiased_exp_het_per_var(
            chunk, pops=of_each_pop, min_num_samples=min_num_individuals
        )["exp_het"]
        for pop in of_each_pop:
            of_the_pop = values[pop].to_numpy()
            of_the_pop = of_the_pop[~numpy.isnan(of_the_pop)]
            on_each_edge[pop] += numpy.isclose(
                of_the_pop[:, numpy.newaxis], edges[numpy.newaxis, :], rtol=0, atol=1e-9
            ).sum(axis=0)
    return on_each_edge


def _compare_the_unbiased_histogram(our_counts, their_counts, on_each_edge) -> None:
    """The histogram of the unbiased expected heterozygosity of one
    population, which agrees with pyNei's but for the variants that lie on
    an edge.

    Each of them moves one count from one of the two bins that share its
    edge to the other, and nothing else does, so the variants below each
    edge can differ by the ones that lie on that edge and by no more, and
    the counts of the two histograms add up to the same number.
    """
    assert our_counts.sum() == their_counts.sum()
    ours_below, theirs_below = numpy.cumsum(our_counts), numpy.cumsum(their_counts)
    # The first edge and the last have no bin on one side, so no variant
    # crosses them: the values below the first and above the last are in no
    # bin in either library.
    for edge in range(1, len(on_each_edge) - 1):
        assert abs(int(ours_below[edge - 1]) - int(theirs_below[edge - 1])) <= int(
            on_each_edge[edge]
        ), (edge, ours_below[edge - 1], theirs_below[edge - 1], on_each_edge[edge])


def _compare_with_pynei(ours, plain, unbiased, pop_names, on_each_edge) -> None:
    """popnei's result against the two runs of pyNei on the same variants.

    `plain` is pyNei's result with `unbiased_exp_het=False`, whose `exp_het`
    is popnei's `exp_het`, and `unbiased` its result with the argument true,
    whose `exp_het` is popnei's `unbiased_exp_het`. The histogram counts and
    the three counts of the polymorphism ratio have to be equal, and the
    means and the two ratios equal within 1e-12 relative.

    The unbiased expected heterozygosity is the one exception, and
    `on_each_edge` is how many variants of each population lie on each edge
    of its histogram, which is what its counts are allowed to differ by.
    """
    theirs_of = {
        "obs_het": plain.obs_het,
        "maf": plain.maf,
        "exp_het": plain.exp_het,
        "unbiased_exp_het": unbiased.exp_het,
    }
    for stat in DISTRIBS:
        ours_of_the_stat = getattr(ours, stat)
        theirs_of_the_stat = theirs_of[stat]
        numpy.testing.assert_array_equal(
            ours_of_the_stat.hist_bin_edges, theirs_of_the_stat.hist_bin_edges
        )
        for pop in pop_names:
            assert _the_same_number(
                ours_of_the_stat.mean[pop], theirs_of_the_stat.mean[pop]
            ), (stat, pop, ours_of_the_stat.mean[pop], theirs_of_the_stat.mean[pop])
            our_counts = numpy.asarray(ours_of_the_stat.hist_counts[pop]).astype(int)
            their_counts = numpy.asarray(theirs_of_the_stat.hist_counts[pop]).astype(
                int
            )
            if stat == "unbiased_exp_het":
                _compare_the_unbiased_histogram(
                    our_counts, their_counts, on_each_edge[pop]
                )
            else:
                numpy.testing.assert_array_equal(
                    our_counts, their_counts, err_msg=f"{stat} of {pop}"
                )
    ours_poly, theirs_poly = ours.poly_vars_ratio, plain.poly_vars_ratio
    for pop in pop_names:
        for count in ("num_poly", "num_variable", "tot_num_variants_with_data"):
            assert int(getattr(ours_poly, count)[pop]) == int(
                getattr(theirs_poly, count)[pop]
            ), (count, pop)
        for ratio in ("poly_ratio", "poly_ratio_over_variables"):
            assert _the_same_number(
                getattr(ours_poly, ratio)[pop], getattr(theirs_poly, ratio)[pop]
            ), (ratio, pop)


def _pyneis_two_runs(path: Path, pops, min_num_individuals: int):
    """pyNei's result on the variants of `path`, with the unbiased expected
    heterozygosity and without it."""
    theirs = vars_from_vcf(path)
    plain = pynei_calc_per_var_distribs(
        theirs, pops=pops, min_num_samples=min_num_individuals, unbiased_exp_het=False
    )
    unbiased = pynei_calc_per_var_distribs(
        theirs, pops=pops, min_num_samples=min_num_individuals, unbiased_exp_het=True
    )
    return plain, unbiased


def test_per_var_distribs_of_the_panel_are_pyneis_over_its_three_populations() -> None:
    """The five statistics of the 1200 variants of the panel, in `p0`, `p1`
    and `p2`, with the default histogram and `min_num_individuals` of 20."""
    ours = calc_per_var_distribs(
        _panel(), pops=PANEL_POPS, min_num_individuals=PANEL_MIN_NUM_INDIVIDUALS
    )
    plain, unbiased = _pyneis_two_runs(PANEL, PANEL_POPS, PANEL_MIN_NUM_INDIVIDUALS)

    assert ours.pass_stats.num_vars == PANEL_NUM_VARS
    on_each_edge = _on_each_bin_edge(
        PANEL, PANEL_POPS, PANEL_MIN_NUM_INDIVIDUALS, ours.maf.hist_bin_edges
    )
    _compare_with_pynei(ours, plain, unbiased, list(PANEL_POPS), on_each_edge)


def test_per_var_distribs_of_many_vcf_are_pyneis_over_its_two_populations() -> None:
    """The same on the 500 variants of `many.vcf`, whose half called
    genotypes and third alleles the panel has none of, with
    `min_num_individuals` of 5."""
    ours = calc_per_var_distribs(
        _many(), pops=MANY_POPS, min_num_individuals=MANY_MIN_NUM_INDIVIDUALS
    )
    plain, unbiased = _pyneis_two_runs(MANY, MANY_POPS, MANY_MIN_NUM_INDIVIDUALS)

    assert ours.pass_stats.num_vars == MANY_NUM_VARS
    on_each_edge = _on_each_bin_edge(
        MANY, MANY_POPS, MANY_MIN_NUM_INDIVIDUALS, ours.maf.hist_bin_edges
    )
    _compare_with_pynei(ours, plain, unbiased, list(MANY_POPS), on_each_edge)


@pytest.mark.parametrize(
    ("path", "min_num_individuals"),
    [
        (PANEL, PANEL_MIN_NUM_INDIVIDUALS),
        (MANY, MANY_MIN_NUM_INDIVIDUALS),
    ],
)
def test_per_var_distribs_with_no_pops_are_pyneis_over_every_individual(
    path: Path, min_num_individuals: int
) -> None:
    """With no `pops` there is one population, named `pop`, of every
    individual, in both libraries."""
    ours = calc_per_var_distribs(
        _panel() if path == PANEL else _many(),
        min_num_individuals=min_num_individuals,
    )
    plain, unbiased = _pyneis_two_runs(path, None, min_num_individuals)

    assert list(ours.maf.mean.index) == ["pop"]
    on_each_edge = _on_each_bin_edge(
        path, None, min_num_individuals, ours.maf.hist_bin_edges
    )
    _compare_with_pynei(ours, plain, unbiased, ["pop"], on_each_edge)


def test_per_var_distribs_leave_a_population_of_15_without_a_value_by_default() -> None:
    """`min_num_individuals` is 20 when nobody gives it, so a population of
    15 diploid individuals, which calls 15 genotypes at most at any variant,
    has no value anywhere: every mean is NaN, every histogram is empty and
    the counts of the polymorphism ratio are 0."""
    fifteen = {"small": PANEL_POPS["p0"][:15]}
    ours = calc_per_var_distribs(_panel(), pops=fifteen)

    for stat in DISTRIBS:
        distrib = getattr(ours, stat)
        assert math.isnan(float(distrib.mean["small"])), stat
        assert int(numpy.asarray(distrib.hist_counts["small"]).sum()) == 0, stat
    poly = ours.poly_vars_ratio
    assert int(poly.num_poly["small"]) == 0
    assert int(poly.num_variable["small"]) == 0
    assert int(poly.tot_num_variants_with_data["small"]) == 0
    assert math.isnan(float(poly.poly_ratio["small"]))
    assert math.isnan(float(poly.poly_ratio_over_variables["small"]))


def test_per_var_distribs_calculates_the_five_statistics_by_default() -> None:
    """Every statistic is there when `stats` is not given, and each one is
    keyed by the names of the populations."""
    ours = calc_per_var_distribs(
        _many(), pops=MANY_POPS, min_num_individuals=MANY_MIN_NUM_INDIVIDUALS
    )

    for stat in DISTRIBS:
        assert isinstance(getattr(ours, stat), StatsDistrib), stat
    assert isinstance(ours.poly_vars_ratio, PolyVarsStats)
    assert list(ours.maf.mean.index) == ["popA", "popB"]


def test_per_var_distribs_give_every_count_as_a_signed_number() -> None:
    """The counts are signed 64 bit integers, as pyNei's are, so that the
    difference of two of them is a negative number and not 1.8e19.

    Over the one population of `many.vcf`, 16 more variants vary than are
    polymorphic at the threshold of 0.95, which pyNei gives as -16 and an
    unsigned subtraction as 18446744073709551600.
    """
    ours = calc_per_var_distribs(_many(), min_num_individuals=MANY_MIN_NUM_INDIVIDUALS)

    poly = ours.poly_vars_ratio
    assert int((poly.num_poly - poly.num_variable)["pop"]) == -16
    for counts in (poly.num_poly, poly.num_variable, poly.tot_num_variants_with_data):
        assert counts.dtype == numpy.int64
    for stat in DISTRIBS:
        counts = getattr(ours, stat).hist_counts["pop"]
        assert counts.dtype == numpy.int64, stat
        assert int((counts - int(counts.max())).min()) < 0, stat


def test_per_var_distribs_calculates_only_the_statistics_asked_for() -> None:
    """A statistic nobody asked for is `None` in the result."""
    ours = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.MAF,),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
    )

    assert isinstance(ours.maf, StatsDistrib)
    assert ours.obs_het is None
    assert ours.exp_het is None
    assert ours.unbiased_exp_het is None
    assert ours.poly_vars_ratio is None


def test_per_var_distribs_takes_one_statistic_on_its_own() -> None:
    """One member of `PerVarStat` written without its tuple asks for that
    statistic, and not for the letters of its name."""
    ours = calc_per_var_distribs(
        _many(), stats=PerVarStat.MAF, min_num_individuals=MANY_MIN_NUM_INDIVIDUALS
    )

    assert isinstance(ours.maf, StatsDistrib)
    assert ours.exp_het is None


def test_per_var_distribs_refuses_a_statistic_written_as_a_string() -> None:
    """`stats` takes the members of `PerVarStat` alone, so that a name with
    a typo in it cannot pass: a string is a `TypeError`, where pyNei takes
    the names as well and refuses an unknown one with a `ValueError`.

    One name written without its comma is refused as the name the user
    wrote, and not as its first letter: a string is a sequence of its
    letters, so `stats="maf"` would otherwise be read as `m`, `a` and `f`.
    """
    with pytest.raises(TypeError, match="'maf', a str, is not one of them"):
        calc_per_var_distribs(_many(), stats="maf")
    with pytest.raises(TypeError, match="'maf', a str, is not one of them"):
        calc_per_var_distribs(_many(), stats=("maf",))
    with pytest.raises(TypeError, match="'mafs', a str, is not one of them"):
        calc_per_var_distribs(_many(), stats=("mafs",))


def test_per_var_distribs_refuse_a_stats_that_names_nothing_at_all() -> None:
    """What is no sequence of members names no statistic, and the message
    says which argument it was and what was given, where Python's own says
    only that an int cannot be iterated over."""
    with pytest.raises(TypeError, match="`stats` is a sequence.*and 5, a int"):
        calc_per_var_distribs(_many(), stats=5)


def test_per_var_distribs_refuses_no_statistic_at_all() -> None:
    """A `stats` with nothing in it asks for a result with nothing in it."""
    with pytest.raises(ValueError, match="statistic"):
        calc_per_var_distribs(_many(), stats=())


def test_per_var_distribs_of_one_pass_are_those_of_one_statistic_at_a_time() -> None:
    """The five statistics of one pass are the five of five passes of one
    each: asking for fewer is a saving of work and changes no value."""
    together = calc_per_var_distribs(
        _many(),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        hist_kwargs={"num_bins": 10},
    )
    for stat in PerVarStat:
        alone = calc_per_var_distribs(
            _many(),
            stats=(stat,),
            pops=MANY_POPS,
            min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
            hist_kwargs={"num_bins": 10},
        )
        one, in_the_pass = getattr(alone, str(stat)), getattr(together, str(stat))
        if stat is PerVarStat.POLY_VARS_RATIO:
            numpy.testing.assert_array_equal(one.num_poly, in_the_pass.num_poly)
            numpy.testing.assert_allclose(one.poly_ratio, in_the_pass.poly_ratio)
        else:
            numpy.testing.assert_allclose(one.mean, in_the_pass.mean)
            numpy.testing.assert_array_equal(
                one.hist_counts.to_numpy(), in_the_pass.hist_counts.to_numpy()
            )


def test_per_var_distribs_arguments_of_one_statistic_change_no_other() -> None:
    """`ploidy` is the two expected heterozygosities' and `poly_threshold`
    the polymorphism ratio's: neither touches the major allele frequency."""
    at_two = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.EXP_HET, PerVarStat.MAF),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
    )
    at_four = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.EXP_HET, PerVarStat.MAF),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        ploidy=4,
    )
    assert not numpy.allclose(at_two.exp_het.mean, at_four.exp_het.mean)
    numpy.testing.assert_allclose(at_two.maf.mean, at_four.maf.mean)

    strict = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.POLY_VARS_RATIO, PerVarStat.MAF),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        poly_threshold=0.51,
    )
    lax = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.POLY_VARS_RATIO, PerVarStat.MAF),
        pops=MANY_POPS,
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        poly_threshold=0.99,
    )
    assert (strict.poly_vars_ratio.num_poly < lax.poly_vars_ratio.num_poly).all()
    numpy.testing.assert_allclose(strict.maf.mean, lax.maf.mean)


def test_per_var_distribs_of_a_pass_whose_filter_kept_no_variant_is_refused() -> None:
    """A mean over no variant is no number, so the pass is an error, and the
    message says that the steps kept none and what each filter counted."""
    variants = _many()
    # No variant with a called allele has a major allele frequency of 0, and
    # one without a called allele is not kept either.
    variants.filter_by_maf(0)

    with pytest.raises(ValueError, match="the pass gave no variant") as refusal:
        calc_per_var_distribs(variants)
    assert "the `maf` filter was given 500 and kept 0" in str(refusal.value)


def test_per_var_distribs_of_a_source_with_no_variant_is_refused(write_vcf) -> None:
    """A VCF with a header and no data line holds no variant, which the
    message says apart from the steps keeping none."""
    with pytest.raises(ValueError, match="its source holds none"):
        calc_per_var_distribs(open_vcf(write_vcf([])))


def test_per_var_distribs_does_not_change_hist_kwargs() -> None:
    """The dict of the histogram is read and not written into."""
    hist_kwargs: dict = {}
    calc_per_var_distribs(
        _many(), min_num_individuals=MANY_MIN_NUM_INDIVIDUALS, hist_kwargs=hist_kwargs
    )
    assert hist_kwargs == {}

    hist_kwargs = {"num_bins": 4}
    calc_per_var_distribs(
        _many(), min_num_individuals=MANY_MIN_NUM_INDIVIDUALS, hist_kwargs=hist_kwargs
    )
    assert hist_kwargs == {"num_bins": 4}


def test_per_var_distribs_takes_the_range_and_the_number_of_bins() -> None:
    """Two bins from 0 to 0.5 have the edges 0, 0.25 and 0.5."""
    ours = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.OBS_HET,),
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        hist_kwargs={"num_bins": 2, "range": (0, 0.5)},
    )

    numpy.testing.assert_array_equal(ours.obs_het.hist_bin_edges, [0, 0.25, 0.5])


def test_per_var_distribs_give_bin_edges_that_nothing_writes_into() -> None:
    """The four distributions of one result share one array of edges, as
    pyNei's do, so a number written into the edges of one statistic would be
    in the edges of the other three: the array is read only."""
    ours = calc_per_var_distribs(_many(), min_num_individuals=MANY_MIN_NUM_INDIVIDUALS)

    assert ours.maf.hist_bin_edges is ours.obs_het.hist_bin_edges
    with pytest.raises(ValueError, match="read-only"):
        ours.maf.hist_bin_edges[0] = 0.5
    assert float(ours.obs_het.hist_bin_edges[0]) == 0


def test_per_var_distribs_makes_logarithmic_bins_span_the_range() -> None:
    """Four bins of equal ratio from 0.01 to 100 have the edges 0.01, 0.1,
    1, 10 and 100."""
    ours = calc_per_var_distribs(
        _many(),
        stats=(PerVarStat.OBS_HET,),
        min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        hist_kwargs={"range": (0.01, 100), "num_bins": 4, "bin_type": "logarithmic"},
    )

    numpy.testing.assert_allclose(
        ours.obs_het.hist_bin_edges, [0.01, 0.1, 1, 10, 100], rtol=1e-12
    )


def test_per_var_distribs_refuse_a_histogram_of_no_bin() -> None:
    """A histogram counts the variants that fall in each of its bins, so one
    with no bin counts nothing: pyNei hands `num_bins` to `linspace`, which
    gives one edge and a histogram that no value falls in."""
    with pytest.raises(ValueError, match="0 bins"):
        calc_per_var_distribs(_many(), hist_kwargs={"num_bins": 0})


def test_per_var_distribs_refuse_more_bins_than_a_machine_counts_in() -> None:
    """A bin is a count of 8 bytes for each population and each statistic,
    once in the pass and once more in every chunk of rows a thread reads, so
    there are 100000 bins at most, far above the tens a person reads. Above
    that bound 2**60 bins were a `PanicException` of a capacity that
    overflowed, which `except Exception` does not catch."""
    for num_bins in (100001, 2**60):
        with pytest.raises(ValueError, match="100000 at most"):
            calc_per_var_distribs(_many(), hist_kwargs={"num_bins": num_bins})


def test_per_var_distribs_refuse_a_range_that_does_not_go_up() -> None:
    """The bins divide the range, so it runs from a number up to a larger
    one: two ends that are equal or the wrong way round give no bin to
    divide, and one that is NaN or infinite leaves every edge NaN."""
    for hist_range in ((1, 1), (1, 0), (0, float("nan")), (0, float("inf"))):
        with pytest.raises(ValueError, match="runs from a number up to a larger one"):
            calc_per_var_distribs(_many(), hist_kwargs={"range": hist_range})


def test_per_var_distribs_refuse_a_range_whose_ends_are_too_far_apart() -> None:
    """The width of a bin is the distance between the two ends over the
    bins, and two ends that are each a number can be further apart than a
    float64 goes: from -1e308 to 1e308 in 4 bins the edges are NaN, infinite,
    infinite, infinite and 1e308, which do not go up, so every value would
    land in the first bin. numpy refuses the same range."""
    with pytest.raises(ValueError, match="above the largest float64"):
        calc_per_var_distribs(_many(), hist_kwargs={"range": (-1e308, 1e308)})


def test_per_var_distribs_refuse_an_unknown_bin_type() -> None:
    """The bins are `linear` or `logarithmic`, and popnei spells the first
    one as English does, where pyNei spells it `lineal`."""
    for bin_type in ("quadratic", "lineal"):
        with pytest.raises(ValueError, match="linear"):
            calc_per_var_distribs(_many(), hist_kwargs={"bin_type": bin_type})


def test_per_var_distribs_refuse_a_logarithmic_range_that_starts_at_zero_or_below() -> (
    None
):
    """Each edge of bins of equal ratio is the one before it times a fixed
    factor, and no factor takes 0 anywhere."""
    for start in (0, -1):
        with pytest.raises(ValueError, match="equal ratio"):
            calc_per_var_distribs(
                _many(),
                hist_kwargs={"range": (start, 10), "bin_type": "logarithmic"},
            )


def test_per_var_distribs_refuse_a_key_of_hist_kwargs_that_is_not_one_of_the_three() -> (
    None
):
    """`range`, `num_bins` and `bin_type` are the histogram, and a key that
    is none of them is a name with a typo in it, which would leave the
    default silently in its place."""
    with pytest.raises(ValueError, match="num_bins"):
        calc_per_var_distribs(_many(), hist_kwargs={"nbins": 4})


@pytest.mark.parametrize("stat", list(PerVarStat))
def test_per_var_distribs_key_every_statistic_by_the_population_names(
    stat: PerVarStat,
) -> None:
    """Each of the five is calculated per population, keyed by the names the
    user gave, and each refuses an individual that is not in the pass."""
    ours = getattr(
        calc_per_var_distribs(
            _many(),
            stats=(stat,),
            pops=MANY_POPS,
            min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
        ),
        str(stat),
    )
    per_pop = ours.num_poly if stat is PerVarStat.POLY_VARS_RATIO else ours.mean
    assert sorted(per_pop.index) == ["popA", "popB"]

    with pytest.raises(ValueError, match="not an individual of the variants"):
        calc_per_var_distribs(_many(), stats=(stat,), pops={"popA": ["nobody"]})


def test_per_var_distribs_keep_the_populations_in_the_order_of_the_keys() -> None:
    """Every statistic holds its values in the order the keys of `pops`
    iterate in, which in Python is the order they were written in."""
    backwards = {"popB": MANY_POPS["popB"], "popA": MANY_POPS["popA"]}
    ours = calc_per_var_distribs(
        _many(), pops=backwards, min_num_individuals=MANY_MIN_NUM_INDIVIDUALS
    )

    for stat in DISTRIBS:
        distrib = getattr(ours, stat)
        assert list(distrib.mean.index) == ["popB", "popA"]
        assert list(distrib.hist_counts.columns) == ["popB", "popA"]
    assert list(ours.poly_vars_ratio.num_poly.index) == ["popB", "popA"]


def test_per_var_distribs_refuse_an_individual_named_twice_in_one_population() -> None:
    """pyNei counts that individual twice; popnei refuses it, since a count
    that is wrong and says nothing is what it never gives."""
    with pytest.raises(ValueError, match="named twice"):
        calc_per_var_distribs(_many(), pops={"popA": ["ind00", "ind01", "ind00"]})


def test_per_var_distribs_refuse_a_population_with_no_individual() -> None:
    """Every statistic of a population is calculated over its individuals,
    so a population with none has no value for any of them."""
    with pytest.raises(ValueError, match="names no individual"):
        calc_per_var_distribs(_many(), pops={"popA": ["ind00"], "empty": []})


def test_per_var_distribs_refuse_pops_with_no_population() -> None:
    """A result holds one value per population, so `pops` with none would
    leave one with nothing in it."""
    with pytest.raises(ValueError, match="names no population"):
        calc_per_var_distribs(_many(), pops={})


def test_per_var_distribs_refuse_one_name_written_as_a_string() -> None:
    """The individuals of a population are a sequence of names, and one name
    written without its comma is the letters of that name."""
    with pytest.raises(ValueError, match="sequence of the names"):
        calc_per_var_distribs(_many(), pops={"popA": "ind00"})


def test_per_var_distribs_refuse_a_min_num_individuals_that_is_no_whole_number() -> (
    None
):
    """It is how many called genotypes a population needs at a variant, so
    it is a whole number and 0 or more, and what is none says so under the
    name of the argument, where pyo3 says only that a float cannot be
    interpreted as an integer.

    A truth value is a whole number in Python, so `True` would ask for one
    called genotype without a word, which is not what whoever wrote it
    meant.
    """
    with pytest.raises(TypeError, match="`min_num_individuals` is 3.1"):
        calc_per_var_distribs(_many(), min_num_individuals=3.1)
    with pytest.raises(TypeError, match="`min_num_individuals` is True"):
        calc_per_var_distribs(_many(), min_num_individuals=True)
    with pytest.raises(ValueError, match="min_num_individuals"):
        calc_per_var_distribs(_many(), min_num_individuals=-1)


def test_per_var_distribs_refuse_a_hist_kwargs_that_is_no_dict() -> None:
    """The histogram is a dict of `range`, `num_bins` and `bin_type`, and
    the message of what is none names the argument and what was given, where
    a list is read key by key and asked for a `get` it has not got and a
    number cannot be iterated over at all."""
    with pytest.raises(TypeError, match="`hist_kwargs` is 5, a int"):
        calc_per_var_distribs(_many(), hist_kwargs=5)
    with pytest.raises(TypeError, match=r"`hist_kwargs` is \['range'\], a list"):
        calc_per_var_distribs(_many(), hist_kwargs=["range"])


def test_per_var_distribs_refuse_a_range_that_is_not_two_ends() -> None:
    """`range` is the two ends of the histogram, and the message of what is
    not two of something names the argument, where pyo3 says only that it
    expected a tuple of length 2."""
    with pytest.raises(TypeError, match=r"`hist_kwargs\['range'\]` is the two ends"):
        calc_per_var_distribs(_many(), hist_kwargs={"range": (0, 1, 2)})
    with pytest.raises(TypeError, match=r"`hist_kwargs\['range'\]` is the two ends"):
        calc_per_var_distribs(_many(), hist_kwargs={"range": 1})


def test_per_var_distribs_refuse_a_ploidy_that_no_statistic_is_raised_to() -> None:
    """`ploidy` is what the allele frequencies of the two expected
    heterozygosities are raised to, 1 at least and 255 at most, the largest
    ploidy a reader gives, and the message names the argument the user
    wrote, where the core names the exponent of a statistic of one
    variant."""
    for ploidy in (0, 256):
        with pytest.raises(ValueError, match=f"`ploidy` is {ploidy}"):
            calc_per_var_distribs(
                _many(),
                min_num_individuals=MANY_MIN_NUM_INDIVIDUALS,
                ploidy=ploidy,
            )


def test_per_var_distribs_refuse_a_poly_threshold_that_is_no_frequency() -> None:
    """A major allele frequency lies between 0 and 1, so a threshold outside
    that makes every variant polymorphic or none."""
    with pytest.raises(TypeError, match="poly_threshold"):
        calc_per_var_distribs(_many(), poly_threshold="0.5")
    for threshold in (-0.1, 1.5, float("nan")):
        with pytest.raises(ValueError, match="poly_threshold"):
            calc_per_var_distribs(_many(), poly_threshold=threshold)


def test_per_var_distribs_give_the_counts_of_the_pass_and_of_its_filters() -> None:
    """The result carries how many variants the pass gave and what each
    filter of the `Variants` was given and kept, in the order of the steps.

    The filter of individuals takes no variant away and has no counts, and a
    missing data filter at 0 over `ind05`, `ind00` and `ind49` keeps 423 of
    the 500 variants of `many.vcf`, which is what the pass then gives.
    """
    ours = calc_per_var_distribs(
        _panel(), pops=PANEL_POPS, min_num_individuals=PANEL_MIN_NUM_INDIVIDUALS
    )
    assert ours.pass_stats.num_vars == PANEL_NUM_VARS
    assert ours.pass_stats.filtering == {}

    variants = _many()
    variants.filter_individuals(("ind05", "ind00", "ind49"))
    variants.filter_by_missing_data(0)
    filtered = calc_per_var_distribs(variants, min_num_individuals=1)

    assert filtered.pass_stats.num_vars == 423
    assert list(filtered.pass_stats.filtering) == ["missing_data"]
    assert filtered.pass_stats.filtering["missing_data"].vars_processed == 500
    assert filtered.pass_stats.filtering["missing_data"].vars_kept == 423


# What a subprocess runs to print the result of one pass over `many.vcf`,
# with the threads of rayon set by `RAYON_NUM_THREADS` in its environment.
# There is no `num_threads` argument, so the pool is what says how many
# threads a pass uses, and it is built once per process.
_THE_PASS_IN_A_PROCESS = """
import json, sys
from popnei import calc_per_var_distribs, open_vcf

path = sys.argv[1]
pops = {
    "popA": [f"ind{idx:02d}" for idx in range(20)],
    "popB": [f"ind{idx:02d}" for idx in range(20, 50)],
}
distribs = calc_per_var_distribs(
    open_vcf(path, only_passed=False), pops=pops, min_num_individuals=5
)
print(
    json.dumps(
        {
            "means": {
                stat: [float(value) for value in getattr(distribs, stat).mean]
                for stat in ("obs_het", "maf", "exp_het", "unbiased_exp_het")
            },
            "hist": {
                stat: [
                    int(count)
                    for count in getattr(distribs, stat).hist_counts.to_numpy().ravel()
                ]
                for stat in ("obs_het", "maf", "exp_het", "unbiased_exp_het")
            },
            "poly": [int(count) for count in distribs.poly_vars_ratio.num_poly],
            "poly_ratio": [
                float(ratio) for ratio in distribs.poly_vars_ratio.poly_ratio
            ],
        }
    )
)
"""


def _the_pass_with_threads(num_threads: int) -> dict:
    """The result of the pass over `many.vcf` in a process whose rayon pool
    has `num_threads` threads."""
    environment = dict(os.environ, RAYON_NUM_THREADS=str(num_threads))
    run = subprocess.run(
        [sys.executable, "-c", _THE_PASS_IN_A_PROCESS, str(MANY)],
        capture_output=True,
        text=True,
        check=True,
        env=environment,
    )
    return json.loads(run.stdout)


def test_per_var_distribs_are_the_same_in_pools_of_one_and_of_four_threads() -> None:
    """How many threads rayon has changes no histogram count and no count of
    the polymorphism ratio, and no mean beyond 1e-12 relative: the counts
    are integers and the sums of the means are added in an order that
    depends on the threads."""
    of_one = _the_pass_with_threads(1)
    of_four = _the_pass_with_threads(4)

    assert of_one["hist"] == of_four["hist"]
    assert of_one["poly"] == of_four["poly"]
    for stat in DISTRIBS:
        for ours, theirs in zip(
            of_one["means"][stat], of_four["means"][stat], strict=True
        ):
            assert _the_same_number(ours, theirs), stat
    for ours, theirs in zip(of_one["poly_ratio"], of_four["poly_ratio"], strict=True):
        assert _the_same_number(ours, theirs)
