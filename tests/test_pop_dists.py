"""The distances between populations, from Python.

`docs/specs/dists.md` has the seven measures, the `PopDists` they come in and
the numbers this file asserts. The literals are the spec's, and the files the
reference programs wrote them into are in `tests/reference/pop_dists/`, which
`tests/reference/pop_dists/make_reference.py` writes.

Two programs are compared with here. plink2 v2.0.0-a.7.7 gives Hudson's F_ST
of both panels, which it prints to six digits, so the comparison is within
1e-6 absolute. ADMIXTOOLS 2.0.10 gives f_2 and its jackknife standard error
of the biallelic panel, which it writes to seventeen digits, so that
comparison is within 1e-12 relative.

ADMIXTOOLS was run on the biallelic panel at three lengths of the resampling
groups and `calc_pop_dists` refuses a pass of fewer than 20 of them, so of
the three runs only the one at 55 000 base pairs, which cuts the panel into
22 groups, can be asked for through this package; the other two, at 100 000
and at 250 000, are cargo tests of `crates/popnei/src/pop_dists.rs`.

Five of the seven measures, the chord distance, Nei's D_A, Jost's D, Nei's
G_ST and the standardized G''_ST, are not calculated yet: the work packages 2
and 3 of `docs/plans/dists-pops.md` add them, and asking for one of them is a
`ValueError` until then, which `test_a_measure_that_is_not_written_yet_is_refused`
holds to.
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

# The two measures that work package 1 of `docs/plans/dists-pops.md`
# calculates, and the five that the work packages 2 and 3 add.
MEASURES_THAT_ARE_WRITTEN = ("fst", "f2")
MEASURES_THAT_ARE_NOT_WRITTEN_YET = ("chord", "da", "dest", "gst", "gst_standardized")


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
    for measure in MEASURES_THAT_ARE_NOT_WRITTEN_YET:
        assert getattr(of_both, measure) is None


def test_a_measure_that_is_not_written_yet_is_refused() -> None:
    """The five measures the work packages 2 and 3 of
    `docs/plans/dists-pops.md` add, which have no value today.

    A user who asks for one gets a `ValueError` that names it and the two
    that are calculated, and not a vector of NaN that says nothing about
    itself. `measures=None`, which asks for all seven, is refused for the
    same reason until those work packages are done.
    """
    for measure in MEASURES_THAT_ARE_NOT_WRITTEN_YET:
        with pytest.raises(ValueError) as refusal:
            calc_pop_dists(
                open_vcf(PANEL),
                PANEL_POPS,
                jackknife_group=None,
                measures=(measure,),
            )
        said = str(refusal.value)
        assert measure in said
        for written in MEASURES_THAT_ARE_WRITTEN:
            assert written in said

    with pytest.raises(ValueError):
        calc_pop_dists(open_vcf(PANEL), PANEL_POPS, jackknife_group=None)


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


def test_a_variants_is_what_the_calculation_takes() -> None:
    """The path of the VCF in the place of the `Variants` is the mistake that
    is easiest to make, and what it gave was the `AttributeError` of an
    object with no source inside it."""
    with pytest.raises(TypeError) as refusal:
        calc_pop_dists(PANEL, PANEL_POPS, jackknife_group=None, measures=("fst",))

    assert "open_vcf" in str(refusal.value)
