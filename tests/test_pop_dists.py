"""The distances between populations, from Python.

`docs/specs/dists.md` has the seven measures, the `PopDists` they come in and
the numbers this file asserts. The literals are the spec's, and the files the
reference programs wrote them into are in `tests/reference/pop_dists/`, which
`tests/reference/pop_dists/make_reference.py` writes.

Four programs are compared with here. plink2 v2.0.0-a.7.7 gives Hudson's
F_ST of both panels, which it prints to six digits, so the comparison is
within 1e-6 absolute. ADMIXTOOLS 2.0.10 gives f_2 and its jackknife standard
error of the biallelic panel, which it writes to seventeen digits, so that
comparison is within 1e-12 relative. pyNei, at the commit `pyproject.toml`
names, gives Jost's D of both panels and is the one program that computes
the estimator popnei computes, so it is run here and its numbers are matched
within 1e-12 relative. mmod 1.3.3 under R 4.6.1 gives Jost's D, Nei's G_ST
and the standardized G''_ST of both panels with another estimator of each,
so those three are an agreement within 5e-4 and not an equality.

ADMIXTOOLS was run on the biallelic panel at three lengths of the resampling
groups and `calc_pop_dists` refuses a pass of fewer than 20 of them, so of
the three runs only the one at 55 000 base pairs, which cuts the panel into
22 groups, can be asked for through this package; the other two, at 100 000
and at 250 000, are cargo tests of `crates/popnei/src/pop_dists.rs`.

All seven measures are calculated, so the refusal of one that popnei has no
value for has nothing left to refuse, which
`test_no_measure_of_the_seven_is_refused` holds to.
"""

import math
from pathlib import Path

import numpy
import pandas
import pytest
from popnei import (
    Distances,
    PopDistMeasure,
    PopDists,
    calc_pairwise_kosman_dists,
    calc_pop_dists,
    open_vcf,
)
from popnei.pop_dists import _MEASURES_THAT_HAVE_A_VALUE, _the_measures
from pynei import vars_from_vcf
from pynei.dists import calc_jost_dest_pop_dists

DISTS_REFERENCE_DIR = Path(__file__).parent / "reference" / "dists"
STATS_REFERENCE_DIR = Path(__file__).parent / "reference" / "stats"
POP_DISTS_REFERENCE_DIR = Path(__file__).parent / "reference" / "pop_dists"

# The biallelic panel of "How it is verified" of the spec: 1200 variants of
# 200 diploid individuals over two chromosomes of 600 each, 3 in 100
# genotypes missing whole, in three populations of 48, 68 and 84.
PANEL = DISTS_REFERENCE_DIR / "panel.vcf.gz"
PANEL_NUM_VARS = 1200

# The multiallelic panel of the same section, written for this item: 120
# microsatellite loci of 90 diploid individuals in three populations of 30,
# six alleles to a locus, 4 in 100 genotypes missing whole. It is the panel
# that shows that the arithmetic does not assume two alleles.
MICRO = POP_DISTS_REFERENCE_DIR / "micro.vcf.gz"
MICRO_NUM_VARS = 120

# The three pairs of both panels, in the order of the distance vector, which
# is the order the populations are named in.
PAIRS = ("p0-p1", "p0-p2", "p1-p2")

# Hudson's F_ST of the three pairs of each panel, which plink2 v2.0.0-a.7.7
# printed with `--fst popcat method=hudson`, from the F_ST item of the spec.
PANEL_FST = (0.104962, 0.102736, 0.109621)
MICRO_FST = (0.0642281, 0.0694106, 0.0699954)
# plink2 prints six digits, so half a unit of the last of them is the
# distance the numbers can be apart, which is what the spec asks for.
FST_TOLERANCE = 1e-6

# The length in base pairs that cuts the biallelic panel into 22 resampling
# groups, ten of 55 variants and one of 50 on each of its two chromosomes.
# It is the only one of the three lengths ADMIXTOOLS was run at that reaches
# the 20 groups `calc_pop_dists` asks for.
PANEL_JACKKNIFE_GROUP = 55000
PANEL_NUM_GROUPS = 22
# How many variants each of those groups holds, in the order they were
# started: the panel has its variants 1000 base pairs apart, so 55 of them
# fall in the first 55 000 of a chromosome and the eleventh group of each
# chromosome holds the 50 that are left.
PANEL_VARS_OF_EACH_GROUP = ([55] * 10 + [50]) * 2

# f_2 and its jackknife standard error of the three pairs of the biallelic
# panel at that length, which ADMIXTOOLS 2.0.10 wrote into
# `tests/reference/pop_dists/panel.f2.min20.tsv`, from the f_2 item of the
# spec and from "The standard errors".
PANEL_F2 = (
    0.041181109098151751,
    0.039890789655075108,
    0.042798563747209556,
)
PANEL_F2_STANDARD_ERRORS = (
    0.0020502481330704485,
    0.0016837006366670754,
    0.0019859616713111257,
)
# Both libraries add the same numbers in a different order, so the last bits
# are what they can differ by, which is the tolerance the spec gives.
F2_TOLERANCE = 1e-12

# The threshold of called genotypes at which the pairs of the biallelic panel
# part, from "How it is verified" of the spec: the variants that count are
# 688 for p0-p1 and p0-p2, whose p0 has 48 individuals, and all 1200 for
# p1-p2. Every variant counts for every pair at the default of 20.
PARTING_MIN_NUM_INDIVIDUALS = 47
PANEL_NUM_VARS_OF_EACH_PAIR = (688, 688, 1200)

# Jost's D of the three pairs of each panel at the default threshold of 20
# called genotypes, which pyNei's `calc_jost_dest_pop_dists` gives and which
# the Jost's D item of the spec prints to ten digits. popnei computes pyNei's
# estimator, so the two agree to the last bits of a double.
PANEL_DEST = (0.0635434630, 0.0612981314, 0.0656705213)
MICRO_DEST = (0.1661307946, 0.1822136338, 0.1819999162)

# The same three of the biallelic panel at the 47 where its pairs part. The
# D of p1-p2 is the one above, since that pair keeps all 1200 variants.
PANEL_DEST_AT_THE_PARTING = (0.0595097904, 0.0612873957, 0.0656705213)

# Nei's G_ST and the standardized G''_ST of the three pairs of each panel at
# that same default, which the G_ST item of the spec prints to ten digits as
# popnei's own numbers: no program outside popnei computes this estimator of
# either, so what they are checked against is mmod below.
PANEL_GST = (0.0554614481, 0.0542279036, 0.0580557529)
PANEL_GST_STANDARDIZED = (0.1619596310, 0.1578689665, 0.1682042511)
MICRO_GST = (0.0331595308, 0.0360169279, 0.0363566716)
MICRO_GST_STANDARDIZED = (0.2196573038, 0.2390740032, 0.2393928219)

# pyNei and popnei add the same numbers in the same order, so the last bits
# are what they can differ by, the tolerance the spec gives for the one
# comparison that pins the estimator of Jost's D.
DEST_PYNEI_TOLERANCE = 1e-12

# The four sets of literals above are the spec's ten digits, so half a unit
# of the tenth is how far a value may be from the one written there.
TEN_DIGITS_TOLERANCE = 5e-11

# Jost's D, Nei's G_ST and the standardized G''_ST of the three pairs of each
# panel as `pairwise_D`, `pairwise_Gst_Nei` and `pairwise_Gst_Hedrick` of mmod
# 1.3.3 under R 4.6.1 give them, which
# `tests/reference/pop_dists/panel.mmod.tsv` and `micro.mmod.tsv` hold.
# `pairwise_Gst_Hedrick` computes the standardized G''_ST of Meirmans and
# Hedrick (2011) and not the G'_ST its name suggests, which "How it is
# verified" of the G_ST item of the spec shows from its source.
PANEL_MMOD = {
    "dest": (0.0634704859, 0.0612312792, 0.0656071139),
    "gst": (0.0553896690, 0.0541614926, 0.0579922831),
    "gst_standardized": (0.1617736271, 0.1576967932, 0.1680418436),
}
MICRO_MMOD = {
    "dest": (0.1662094246, 0.1819580763, 0.1816462134),
    "gst": (0.0331789783, 0.0359529575, 0.0362672006),
    "gst_standardized": (0.2197612680, 0.2387386980, 0.2389275805),
}

# mmod computes another estimator of the same three quantities: it leaves the
# observed heterozygosity term out of both corrections and uses 2n/(2n - 1)
# where popnei, which is pyNei and Nei and Chesser (1983), uses n/(n - 1) and
# subtracts H_obs/(2n). So the check is an agreement and not an equality, and
# 5e-4 is the tolerance the two items of the spec give: the furthest of these
# eighteen numbers, a G''_ST of the multiallelic panel, is 4.7e-4 away.
MMOD_TOLERANCE = 5e-4

# The threshold at which the two pairs of p0 count no variant at all: p0 has
# 48 individuals, so it never has 50 called genotypes at a variant, and p1
# and p2, of 68 and 84, have every variant of the panel at it.
NO_VARIANT_MIN_NUM_INDIVIDUALS = 50

# The seven measures, in the order `PopDistMeasure` has them, all of which
# popnei calculates.
EVERY_MEASURE = (
    "fst",
    "f2",
    "chord",
    "da",
    "dest",
    "gst",
    "gst_standardized",
)


def _pops_of(path: Path) -> dict[str, list[str]]:
    """The populations of a file of `IID` and `popcat` columns, which is what
    plink2 reads them from, under their names in order.

    popnei keeps the order of the `pops` dict, so the pairs come out in the
    order of the names here, which is the order the reference programs print
    their pairs in. The individuals of the biallelic panel are interleaved in
    its file, p0 first, then p2, then p1, so reading them in the order they
    appear would give the pairs in another order than the files below.
    """
    pops: dict[str, list[str]] = {}
    for line in path.read_text().splitlines()[1:]:
        name, pop = line.split("\t")
        pops.setdefault(pop, []).append(name)
    return {pop: pops[pop] for pop in sorted(pops)}


PANEL_POPS = _pops_of(STATS_REFERENCE_DIR / "panel_pops.txt")
MICRO_POPS = _pops_of(POP_DISTS_REFERENCE_DIR / "micro_pops.txt")


def _assert_within(
    found, expected, tolerance: float, relative: bool, what: str
) -> None:
    """Every value of `found` within `tolerance` of the one beside it in
    `expected`, relative to the expected value or absolute."""
    assert len(found) == len(expected), what
    for pair, (ours, theirs) in enumerate(zip(found, expected, strict=True)):
        apart = abs(ours - theirs)
        if relative:
            apart = apart / abs(theirs)
        assert apart <= tolerance, f"{what} of the pair {PAIRS[pair]}: {ours} {theirs}"


def test_the_fst_of_both_panels_is_the_one_plink2_gives() -> None:
    """Hudson's F_ST of the three pairs of the biallelic and of the
    multiallelic panel.

    The multiallelic one is the one that says that every allele counts as
    itself: the same formula on the major allele against the rest gives
    0.0707950 for p0 and p1 there, and on the reference allele against the
    rest 0.0598407, neither of which is what plink2 printed.

    The two panels are cut into resampling groups in different ways, a length
    in base pairs for the one that has positions along two chromosomes and
    one group for each locus for the microsatellites, which have no linkage
    to speak of. Neither changes an F_ST: the groups are what the standard
    errors are built from.
    """
    for path, pops, expected, group in (
        (PANEL, PANEL_POPS, PANEL_FST, PANEL_JACKKNIFE_GROUP),
        (MICRO, MICRO_POPS, MICRO_FST, "variant"),
    ):
        dists = calc_pop_dists(
            open_vcf(path), pops, jackknife_group=group, measures=("fst",)
        )

        of_the_panel = f"the F_ST of {path.name}"
        _assert_within(
            dists.fst.dist_vector, expected, FST_TOLERANCE, False, of_the_panel
        )
        assert dists.pops == ("p0", "p1", "p2"), of_the_panel
        assert dists.fst.names == ("p0", "p1", "p2"), of_the_panel
        assert dists.pass_stats.num_vars == (
            PANEL_NUM_VARS if path == PANEL else MICRO_NUM_VARS
        ), of_the_panel


def test_the_f2_of_the_panel_and_its_standard_error_are_admixtools() -> None:
    """f_2 of the three pairs of the biallelic panel and the standard error
    of each, over the 22 groups that a length of 55 000 base pairs cuts the
    panel into.

    It is the one of the three runs of ADMIXTOOLS that is above the 20 groups
    the function asks for, and the standard errors are what checks the
    delete-m jackknife against a program that computes the same one.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=("f2", "fst"),
    )

    _assert_within(
        dists.f2.dist_vector, PANEL_F2, F2_TOLERANCE, True, "the f_2 of the panel"
    )
    _assert_within(
        dists.f2.standard_errors,
        PANEL_F2_STANDARD_ERRORS,
        F2_TOLERANCE,
        True,
        "the standard error of the f_2 of the panel",
    )
    # The F_ST of the same call, which shares the pass and the groups, has a
    # standard error of the same jackknife, and no program prints one to
    # compare it with: what is asserted is that it is there and is a number
    # of the size of an F_ST of these populations.
    assert dists.fst.standard_errors is not None
    assert all(0 < error < 0.01 for error in dists.fst.standard_errors)


def test_the_dest_of_both_panels_is_the_one_pynei_gives() -> None:
    """Jost's D of the three pairs of each panel, against pyNei's
    `calc_jost_dest_pop_dists` run on the same file.

    pyNei is the one program that computes the estimator popnei computes,
    the Nei and Chesser correction that GenAlEx prints, so this is the
    comparison that pins it: the two add the same numbers in the same order
    and agree to the last bits of a double. mmod's D, which the test below
    compares with, is another estimator of the same quantity and is 7.3e-5
    away on the biallelic panel and 3.5e-4 on the multiallelic one, so it
    says that popnei computes Jost's D and not which estimator of it.

    Both panels are read at the default of 20 called genotypes, where every
    variant counts for every pair of both of them. The multiallelic one is
    the one that says that the arithmetic does not assume two alleles: its
    120 loci have six alleles each.
    """
    for path, pops, expected in (
        (PANEL, PANEL_POPS, PANEL_DEST),
        (MICRO, MICRO_POPS, MICRO_DEST),
    ):
        ours = calc_pop_dists(
            open_vcf(path), pops, jackknife_group=None, measures=("dest",)
        )
        theirs = calc_jost_dest_pop_dists(vars_from_vcf(path), pops)

        of_the_panel = f"Jost's D of {path.name}"
        _assert_within(
            ours.dest.dist_vector,
            theirs.dist_vector,
            DEST_PYNEI_TOLERANCE,
            True,
            of_the_panel,
        )
        _assert_within(
            ours.dest.dist_vector,
            expected,
            TEN_DIGITS_TOLERANCE,
            False,
            of_the_panel,
        )


def test_the_dest_where_the_pairs_part_is_over_each_pairs_own_variants() -> None:
    """Jost's D of the biallelic panel at 47 called genotypes, where the two
    pairs of p0 lose the variants at which p0, of 48 individuals, has fewer
    than 47 called.

    Each pair is a mean over its own variants: 688 for p0-p1 and p0-p2 and
    all 1200 for p1-p2, whose D is therefore the one at the default of 20
    and does not move. pyNei drops the same variants by another route, a
    frequency it sets to NaN and a second test on the called genotypes, and
    gives the same three numbers.
    """
    ours = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=None,
        measures=("dest",),
        min_num_individuals=PARTING_MIN_NUM_INDIVIDUALS,
    )
    theirs = calc_jost_dest_pop_dists(
        vars_from_vcf(PANEL),
        PANEL_POPS,
        min_num_samples=PARTING_MIN_NUM_INDIVIDUALS,
    )

    at_the_parting = "Jost's D of the panel at 47 called genotypes"
    _assert_within(
        ours.dest.dist_vector,
        theirs.dist_vector,
        DEST_PYNEI_TOLERANCE,
        True,
        at_the_parting,
    )
    _assert_within(
        ours.dest.dist_vector,
        PANEL_DEST_AT_THE_PARTING,
        TEN_DIGITS_TOLERANCE,
        False,
        at_the_parting,
    )
    assert tuple(ours.num_vars) == PANEL_NUM_VARS_OF_EACH_PAIR


def test_the_dest_the_gst_and_the_gst_standardized_of_both_panels_agree_with_mmod() -> (
    None
):
    """The three measures that come out of the corrected H_S and H_T, of the
    three pairs of both panels, against `pairwise_D`, `pairwise_Gst_Nei` and
    `pairwise_Gst_Hedrick` of mmod 1.3.3.

    mmod computes another estimator of each of the three, so this is an
    agreement within 5e-4 and not an equality, and it says that popnei
    computes these three quantities and not other statistics. What pins the
    estimator is the comparison with pyNei above, exact to 1e-12 relative.
    Tightening this tolerance, or moving the arithmetic towards mmod, breaks
    that one.

    The same call also gives the ten digits the spec prints for popnei's own
    G_ST and G''_ST, which no program outside popnei computes and which are
    here so that a change to the arithmetic that stays inside 5e-4 of mmod
    is still caught.
    """
    for path, pops, of_mmod, of_popnei in (
        (PANEL, PANEL_POPS, PANEL_MMOD, (PANEL_GST, PANEL_GST_STANDARDIZED)),
        (MICRO, MICRO_POPS, MICRO_MMOD, (MICRO_GST, MICRO_GST_STANDARDIZED)),
    ):
        dists = calc_pop_dists(
            open_vcf(path),
            pops,
            jackknife_group=None,
            measures=("dest", "gst", "gst_standardized"),
        )

        for measure, expected in of_mmod.items():
            _assert_within(
                getattr(dists, measure).dist_vector,
                expected,
                MMOD_TOLERANCE,
                False,
                f"the {measure} of {path.name} against mmod",
            )
        for measure, expected in zip(
            ("gst", "gst_standardized"), of_popnei, strict=True
        ):
            _assert_within(
                getattr(dists, measure).dist_vector,
                expected,
                TEN_DIGITS_TOLERANCE,
                False,
                f"the {measure} of {path.name}",
            )


def test_the_groups_the_variants_fell_into_carry_their_chromosome_and_positions() -> (
    None
):
    """The 22 groups of the biallelic panel: 11 on each of its two
    chromosomes, the first of each holding the positions 1000 to 55000, which
    are 55 variants of the 1000 apart the panel has them at.

    A group is anchored on its own first variant, so the first group of chr2
    starts at that chromosome's first position and not where the last group
    of chr1 ended.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=("f2",),
    )

    assert len(dists.group_ids) == PANEL_NUM_GROUPS
    assert dists.group_ids[0] == ("chr1", 1000, 55000)
    assert dists.group_ids[10] == ("chr1", 551000, 600000)
    assert dists.group_ids[11] == ("chr2", 1000, 55000)
    assert dists.group_ids[21] == ("chr2", 551000, 600000)


def test_the_f2_of_each_group_is_a_row_of_f2_groups() -> None:
    """`f2_groups`, the f_2 of every pair within every group, which f_3 and
    f_4 are built from later without reading the genotypes again.

    It is 22 groups by 3 pairs. The f_2 of a pair over the whole panel is the
    mean of its variants, so the f_2 of each group weighted by the variants
    of that group, 55 in twenty of them and 50 in the other two, comes back
    to it: that is what ties the table to the numbers ADMIXTOOLS gave.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=("f2",),
    )

    assert dists.f2_groups.shape == (PANEL_NUM_GROUPS, 3)
    assert dists.f2_groups.dtype == numpy.float64
    of_each_group = numpy.array(PANEL_VARS_OF_EACH_GROUP, dtype=numpy.float64)
    for pair, over_the_panel in enumerate(PANEL_F2):
        of_the_groups = dists.f2_groups[:, pair]
        assert numpy.isfinite(of_the_groups).all(), PAIRS[pair]
        weighted = (of_the_groups * of_each_group).sum() / PANEL_NUM_VARS
        assert abs(weighted - over_the_panel) / over_the_panel < F2_TOLERANCE, PAIRS[
            pair
        ]


def test_no_jackknife_group_gives_no_standard_error_and_no_groups() -> None:
    """`jackknife_group=None`, which asks for no standard errors.

    The measures are the same numbers as with the groups, since the groups
    are only what the standard errors are resampled over, and the result
    carries no `standard_errors`, no `f2_groups` and no group. That the pass
    then asks its reader for no positions is what the cargo test
    `the_pass_asks_for_the_positions_only_where_the_groups_need_them` of
    `crates/popnei/src/pop_dists.rs` holds to: nothing of it is visible from
    Python.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=("f2", "fst")
    )

    _assert_within(
        dists.f2.dist_vector, PANEL_F2, F2_TOLERANCE, True, "the f_2 of the panel"
    )
    assert dists.f2.standard_errors is None
    assert dists.fst.standard_errors is None
    assert dists.f2_groups is None
    assert dists.group_ids == ()


def test_a_call_with_no_jackknife_group_is_refused_by_python_itself() -> None:
    """The length of the resampling groups has no default and the call fails
    without it.

    The right length depends on the linkage disequilibrium of the populations
    being compared, which popnei cannot know, so a user who has not thought
    about it gets no standard error rather than one that looks like the
    others. Python raises for the missing argument itself, and the message
    names it.
    """
    with pytest.raises(TypeError) as refusal:
        calc_pop_dists(open_vcf(PANEL), PANEL_POPS, measures=("fst",))

    assert "jackknife_group" in str(refusal.value)


def test_groups_of_no_base_pairs_are_refused() -> None:
    """A length of 0 base pairs, which is no stretch of a chromosome.

    A user who wants each variant in a group of its own writes `"variant"`
    and one who wants no standard error writes `None`; neither is a length.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL), PANEL_POPS, jackknife_group=0, measures=("fst",)
        )

    said = str(refusal.value)
    assert "`jackknife_group`" in said
    assert "variant" in said


def test_a_jackknife_group_that_is_no_length_and_no_word_is_refused() -> None:
    """The four values beside 0 that are none of the three kinds the argument
    takes, a length in base pairs of 1 or more, `"variant"` and `None`.

    A negative length and the word `"variants"` are the wrong value of a kind
    the argument takes, which is a `ValueError`; 2.5 and `True` are of a kind
    it does not take, which is a `TypeError`. `True` is a whole number in
    Python and would otherwise be a group of one base pair, which says
    nothing about the linkage disequilibrium of the populations being
    compared.
    """
    for given, refused_with in (
        (-1, ValueError),
        ("variants", ValueError),
        (2.5, TypeError),
        (True, TypeError),
    ):
        with pytest.raises(refused_with) as refusal:
            calc_pop_dists(
                open_vcf(PANEL),
                PANEL_POPS,
                jackknife_group=given,
                measures=("fst",),
            )

        said = str(refusal.value)
        assert "`jackknife_group`" in said, given
        assert "variant" in said, given


def test_fewer_than_twenty_groups_are_refused_with_how_many_there_were() -> None:
    """A length of 100 000 base pairs on the biallelic panel, which cuts it
    into 12 groups.

    Each group is left out in turn, so a standard error built from a handful
    of them says more about where the cuts fell than about the populations. A
    user who chose a length too long for their data is told how many groups
    it gave, which is what says how much shorter to make them.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL), PANEL_POPS, jackknife_group=100000, measures=("fst",)
        )

    said = str(refusal.value)
    assert "12" in said
    assert "20" in said


def test_a_source_whose_variants_go_back_is_refused(write_vcf) -> None:
    """Three variants of one chromosome at 1000, 3000 and 2000.

    The groups the standard errors are resampled over are stretches of one
    chromosome, cut by comparing the position of a variant with the first
    position of the group being filled, so a variant that goes back joins
    that group instead of starting one and the groups are not the stretches
    the user asked for. The message names the chromosome and the two
    positions, and the file they were read from.
    """
    path = write_vcf(
        [
            "chr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t3000\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t2000\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
        ]
    )

    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(path),
            {"a": ["ind1"], "b": ["ind2", "ind3"]},
            jackknife_group=1000,
            measures=("fst",),
            min_num_individuals=1,
        )

    said = str(refusal.value)
    assert "chr1" in said
    assert "2000" in said
    assert "3000" in said
    assert str(path) in said


def test_fewer_than_two_populations_are_refused() -> None:
    """One population, which makes no pair.

    Every one of the seven measures is of a pair of populations, so there is
    nothing to give. What the user wrote is wrong whatever file is read, so
    the message names no file.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL),
            {"p0": PANEL_POPS["p0"]},
            jackknife_group=None,
            measures=("fst",),
        )

    said = str(refusal.value)
    assert "pops" in said
    assert str(PANEL) not in said


def test_pops_with_no_population_at_all_is_refused_with_the_rule_of_the_call() -> None:
    """`pops` with nothing in it.

    `pops` is a required argument of `calc_pop_dists`, so the advice the
    statistics give for an empty one, to leave it out for one population of
    every individual, cannot be followed here. What the user has to do is
    name two populations, which is what the message says.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL),
            {},
            jackknife_group=None,
            measures=("fst",),
        )

    said = str(refusal.value)
    assert "pops" in said
    assert "two" in said
    assert "leave" not in said
    assert str(PANEL) not in said


def test_a_name_that_is_not_an_individual_of_the_pass_is_refused() -> None:
    """A population that names an individual the source does not hold.

    The names are looked up among the individuals the pass gives, which are
    those of the source after a filter of individuals when the `Variants` has
    one, so only the pass knows them.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL),
            {"p0": ["s000", "nobody"], "p1": PANEL_POPS["p1"]},
            jackknife_group=None,
            measures=("fst",),
        )

    assert "nobody" in str(refusal.value)


def test_a_pass_that_kept_no_variant_is_refused_with_the_counts_of_its_filters(
    write_vcf,
) -> None:
    """Four variants, each with one genotype of the three missing, and a
    filter of missing data that keeps the variants with no missing genotype.

    The counts of a pass that could not finish are otherwise lost, so the
    message carries them: the source gave four variants and the filter kept
    none of them. A measure over no variant is no number, so it is a wrong
    input and not a vector of NaN.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t./.",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t0/1\t1/1",
            "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
        ]
    )
    variants = open_vcf(path)
    variants.filter_by_missing_data(0)

    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            variants,
            {"one": ["ind1"], "two": ["ind2", "ind3"]},
            jackknife_group=None,
            measures=("fst",),
            min_num_individuals=1,
        )

    assert str(refusal.value) == (
        f"{path}: the pass gave no variant: its source gave 4 and the steps "
        f"kept none of them, the `missing_data` filter was given 4 and kept "
        f"0; a statistic of a pass is calculated over the variants it gives"
    )


def test_the_variants_that_count_are_counted_for_each_pair_on_its_own() -> None:
    """The biallelic panel at a `min_num_individuals` of 47, which is above
    the 42 called genotypes that p0, of 48 individuals, has at its emptiest
    variant and below the 60 and 76 of p1 and p2.

    A variant counts for a pair when both of its populations have at least
    that many called genotypes there, so a variant p0 is short of is lost by
    the two pairs p0 is in and kept by p1-p2: 688, 688 and 1200 of the 1200
    variants. The measures of two pairs are therefore means over different
    variants, which is why each pair carries its own count.

    At the default of 20 every variant counts for every pair, which is the
    other half of the test: a fixture in which the counts are always the same
    number cannot tell a count per pair from one count for all of them.
    """
    parted = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=None,
        measures=("fst",),
        min_num_individuals=PARTING_MIN_NUM_INDIVIDUALS,
    )
    assert tuple(parted.num_vars) == PANEL_NUM_VARS_OF_EACH_PAIR
    assert parted.num_vars.dtype == numpy.int64

    at_the_default = calc_pop_dists(
        open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=("fst",)
    )
    assert tuple(at_the_default.num_vars) == (PANEL_NUM_VARS,) * 3
    # The F_ST of p1-p2 is over the same 1200 variants at both thresholds and
    # is the same number; the two pairs that lost variants moved.
    assert at_the_default.fst.dist_vector[2] == parted.fst.dist_vector[2]
    assert at_the_default.fst.dist_vector[0] != parted.fst.dist_vector[0]


def test_a_pair_with_no_variant_is_nan_and_the_pass_is_not_an_error() -> None:
    """The biallelic panel at a `min_num_individuals` of 50, which p0 never
    reaches: it has 48 individuals, so its largest count of called genotypes
    at a variant is 48.

    The two pairs p0 is in count no variant at all, and the spec asks for no
    value for any measure there, NaN in the distance vector and in the
    standard error, and a `num_vars` of 0. The pass is not an error: p1-p2
    counted every variant and has the f_2 and the standard error that
    ADMIXTOOLS gives for it, which is what says that the two NaN are the two
    pairs and not the whole pass.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=EVERY_MEASURE,
        min_num_individuals=NO_VARIANT_MIN_NUM_INDIVIDUALS,
    )

    assert tuple(dists.num_vars) == (0, 0, PANEL_NUM_VARS)
    for measure in EVERY_MEASURE:
        values = getattr(dists, measure)
        assert math.isnan(values.dist_vector[0]), measure
        assert math.isnan(values.dist_vector[1]), measure
        assert math.isnan(values.standard_errors[0]), measure
        assert math.isnan(values.standard_errors[1]), measure
        assert not math.isnan(values.dist_vector[2]), measure
    _assert_within(
        dists.f2.dist_vector[2:], PANEL_F2[2:], F2_TOLERANCE, True, "the f_2 of p1-p2"
    )
    _assert_within(
        dists.f2.standard_errors[2:],
        PANEL_F2_STANDARD_ERRORS[2:],
        F2_TOLERANCE,
        True,
        "the standard error of the f_2 of p1-p2",
    )
    # The f_2 of the pairs with no variant is NaN within every group too: no
    # variant of theirs fell in any of them.
    assert numpy.isnan(dists.f2_groups[:, :2]).all()
    assert not numpy.isnan(dists.f2_groups[:, 2]).any()
    assert dists.pass_stats.num_vars == PANEL_NUM_VARS


def test_the_pairs_are_in_the_order_the_populations_were_named_in() -> None:
    """The populations stay in the order of the `pops` dict, which is where
    popnei parts from pyNei: pyNei sorts their names.

    No value changes with the order, so the same three numbers come out
    against the other pairs. Named p2, p1, p0, the pairs are p2-p1, p2-p0 and
    p1-p0, and their F_ST is that of p1-p2, of p0-p2 and of p0-p1.
    """
    turned_around = {pop: PANEL_POPS[pop] for pop in ("p2", "p1", "p0")}

    dists = calc_pop_dists(
        open_vcf(PANEL), turned_around, jackknife_group=None, measures=("fst",)
    )

    assert dists.pops == ("p2", "p1", "p0")
    _assert_within(
        dists.fst.dist_vector,
        (PANEL_FST[2], PANEL_FST[1], PANEL_FST[0]),
        FST_TOLERANCE,
        False,
        "the F_ST of the panel named the other way round",
    )


def test_the_counts_of_the_pass_hold_the_variants_and_the_filters() -> None:
    """`pass_stats`, the counts every consumer of a `Variants` gives back:
    how many variants the pass took after the steps, and what each filter was
    given and kept.

    The filter here keeps the variants whose major allele frequency is below
    0.95, and the count of the variants of the pass is what it kept.
    """
    variants = open_vcf(PANEL)
    variants.filter_by_maf(0.95)

    dists = calc_pop_dists(
        variants, PANEL_POPS, jackknife_group=None, measures=("fst",)
    )

    kept = dists.pass_stats.filtering["maf"]
    assert kept.vars_processed == PANEL_NUM_VARS
    assert 0 < kept.vars_kept < PANEL_NUM_VARS
    assert dists.pass_stats.num_vars == kept.vars_kept


def test_each_measure_carries_the_counts_of_the_pass() -> None:
    """The `Distances` of every measure holds the `pass_stats` the `PopDists`
    holds: the measures come out of one pass, so a user who takes one of them
    out of the result keeps the counts of the pass that gave it, as they do
    with the Kosman distances.
    """
    variants = open_vcf(PANEL)
    variants.filter_by_maf(0.95)

    dists = calc_pop_dists(
        variants, PANEL_POPS, jackknife_group=None, measures=("fst", "f2")
    )

    assert dists.fst.pass_stats == dists.pass_stats
    assert dists.f2.pass_stats == dists.pass_stats
    assert dists.pass_stats.filtering["maf"].vars_processed == PANEL_NUM_VARS


def test_the_variants_are_as_they_were_after_the_call() -> None:
    """The call is a consumer of the `Variants`: it makes one pass over the
    source through the steps the `Variants` has, and the `Variants` is as it
    was afterwards, so a second call gives the same numbers."""
    variants = open_vcf(PANEL)

    first = calc_pop_dists(
        variants, PANEL_POPS, jackknife_group=None, measures=("fst",)
    )
    assert variants.steps == ()
    second = calc_pop_dists(
        variants, PANEL_POPS, jackknife_group=None, measures=("fst",)
    )

    assert numpy.array_equal(first.fst.dist_vector, second.fst.dist_vector)


def test_a_measure_that_was_not_asked_for_is_none_in_the_result() -> None:
    """`PopDists` holds one `Distances` for each measure that was asked for,
    under the name of the measure, and `None` for the ones that were not.

    Asking for one measure gives the same number as asking for both: the pass
    is what costs and each measure is a division at the end of it.
    """
    of_one = calc_pop_dists(
        open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=("f2",)
    )
    of_both = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=None,
        measures=(PopDistMeasure.FST, PopDistMeasure.F2),
    )

    assert of_one.fst is None
    assert isinstance(of_one.f2, Distances)
    assert isinstance(of_both.fst, Distances)
    assert numpy.array_equal(of_one.f2.dist_vector, of_both.f2.dist_vector)
    for measure in EVERY_MEASURE:
        if measure not in ("fst", "f2"):
            assert getattr(of_both, measure) is None, measure


def test_one_measure_as_a_bare_string_and_a_repeated_one_are_taken_once() -> None:
    """A member of a `StrEnum` is a string, so one written on its own is a
    sequence of its letters and would be read as six measures of one letter
    each: it is taken as that one measure, as a name written as a string is.

    A name written twice is taken once, so that the result holds one
    `Distances` for it and the pass calculates it once.
    """
    of_a_bare_string = calc_pop_dists(
        open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures="fst"
    )
    of_a_member = calc_pop_dists(
        open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=PopDistMeasure.FST
    )
    of_a_name_twice = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=None,
        measures=("fst", "f2", "fst"),
    )

    assert of_a_bare_string.f2 is None
    assert of_a_member.f2 is None
    assert numpy.array_equal(
        of_a_bare_string.fst.dist_vector, of_a_member.fst.dist_vector
    )
    assert numpy.array_equal(
        of_a_name_twice.fst.dist_vector, of_a_bare_string.fst.dist_vector
    )
    assert of_a_name_twice.f2 is not None
    # A name written twice gives the same result either way, since the second
    # one would overwrite the field of the first with the same numbers: what
    # the result cannot show is that the measure was named once for the pass,
    # so the names it was asked for are read here.
    assert _the_measures(("fst", "f2", "fst")) == ["fst", "f2"]
    assert _the_measures("fst") == ["fst"]


def test_no_measure_of_the_seven_is_refused() -> None:
    """Every member of `PopDistMeasure` has a value, so `measures=None`,
    which asks for all seven, is taken and the refusal of a measure popnei
    has no value for has nothing left to refuse.

    The list of the measures that have one is the core's, and a measure
    added to `PopDistMeasure` without a formula beside it would be refused
    here rather than handed to a user as a vector of NaN.
    """
    assert _MEASURES_THAT_HAVE_A_VALUE == EVERY_MEASURE
    assert _the_measures(None) == list(EVERY_MEASURE)


def test_a_measure_that_is_of_none_of_the_seven_is_refused() -> None:
    """A name written as a string where the members of `PopDistMeasure` go.

    A member of a `StrEnum` is a string, so `"fst"` is taken as the member it
    is; a name with a typo in it is not one of them and is refused with the
    seven, and anything that is no name at all is refused by its type.
    """
    with pytest.raises(ValueError) as refusal:
        calc_pop_dists(
            open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=("fsts",)
        )
    assert "fst" in str(refusal.value)

    with pytest.raises(TypeError):
        calc_pop_dists(open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=(2,))

    with pytest.raises(ValueError):
        calc_pop_dists(open_vcf(PANEL), PANEL_POPS, jackknife_group=None, measures=())


def test_the_measures_come_in_a_frozen_result_of_read_only_arrays() -> None:
    """`PopDists` and the `Distances` in it are frozen dataclasses whose
    arrays cannot be written into: a number written into one of them would be
    a result nobody could tell from a measured one."""
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=("f2",),
    )

    assert isinstance(dists, PopDists)
    with pytest.raises(AttributeError):
        dists.pops = ()
    for array in (
        dists.num_vars,
        dists.f2_groups,
        dists.f2.dist_vector,
        dists.f2.standard_errors,
    ):
        with pytest.raises(ValueError):
            array[0] = 0


def test_the_square_standard_errors_are_the_frame_of_the_vector() -> None:
    """`square_standard_errors`, the N x N frame of the standard errors,
    indexed by the names of the populations on both sides with NaN on its
    diagonal, which is the form a table for a paper takes.

    The Kosman distances have no standard errors, and the method gives `None`
    for them as the field does.
    """
    dists = calc_pop_dists(
        open_vcf(PANEL),
        PANEL_POPS,
        jackknife_group=PANEL_JACKKNIFE_GROUP,
        measures=("f2",),
    )

    square = dists.f2.square_standard_errors()
    assert isinstance(square, pandas.DataFrame)
    assert list(square.index) == ["p0", "p1", "p2"]
    assert list(square.columns) == ["p0", "p1", "p2"]
    assert square.loc["p0", "p1"] == pytest.approx(
        PANEL_F2_STANDARD_ERRORS[0], rel=F2_TOLERANCE
    )
    assert square.loc["p1", "p0"] == square.loc["p0", "p1"]
    assert all(math.isnan(square.iloc[place, place]) for place in range(3))

    assert Distances(dist_vector=[1.0]).square_standard_errors() is None


def test_the_kosman_distances_have_no_standard_errors() -> None:
    """`Distances.standard_errors` is the field the distances between
    populations added, and the Kosman distances between individuals leave it
    `None`: nothing of theirs changes.
    """
    dists = calc_pairwise_kosman_dists(open_vcf(MICRO))

    assert dists.standard_errors is None
    assert dists.square_standard_errors() is None


def test_a_distances_is_built_with_standard_errors_of_its_own() -> None:
    """A `Distances` a user builds themselves from numbers that were
    calculated elsewhere takes the standard errors beside them, one for each
    pair, and a vector of another length is refused."""
    built = Distances(
        dist_vector=[0.1, 0.2, 0.3],
        names=("a", "b", "c"),
        standard_errors=[0.01, 0.02, 0.03],
    )

    assert built.standard_errors.dtype == numpy.float64
    assert not built.standard_errors.flags.writeable
    assert built.square_standard_errors().loc["a", "c"] == 0.02

    with pytest.raises(ValueError) as refusal:
        Distances(dist_vector=[0.1, 0.2, 0.3], standard_errors=[0.01])
    assert "standard_errors" in str(refusal.value)


def test_standard_errors_that_are_not_one_row_of_numbers_are_refused() -> None:
    """A column of standard errors beside a vector of distances.

    The values are as many as the pairs and their shape is not the one a
    result holds, so the message about how many there are, which is the
    other thing that can be wrong with them, would say that three values
    were given for three distances and ask for one for each of them.
    `dist_vector` has a message of its own for the same mistake.
    """
    with pytest.raises(ValueError) as refusal:
        Distances(
            dist_vector=[0.1, 0.2, 0.3],
            standard_errors=[[1.0], [2.0], [3.0]],
        )

    said = str(refusal.value)
    assert "standard_errors" in said
    assert "2 dimensions" in said


def test_a_variants_is_what_the_calculation_takes() -> None:
    """The path of the VCF in the place of the `Variants` is the mistake that
    is easiest to make, and what it gave was the `AttributeError` of an
    object with no source inside it."""
    with pytest.raises(TypeError) as refusal:
        calc_pop_dists(PANEL, PANEL_POPS, jackknife_group=None, measures=("fst",))

    assert "open_vcf" in str(refusal.value)
