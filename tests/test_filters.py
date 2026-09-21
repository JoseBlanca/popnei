"""The three filters from Python: which variants they keep, what they count
and what a `Variants` carries once they are put on it.

`docs/specs/filters.md` has the three filters, the counts of each and the
`Step` that a filter is in a `Variants`. The comparison it asks for is with
pyNei at commit ef0ca6e, which `pyproject.toml` names: `many.vcf` of
`docs/specs/io_vcf.md`, 500 variants of 50 diploid individuals, is read by
both libraries, popnei with `only_passed=False` because pyNei gives every
variant whatever its FILTER says, and the genotypes and the positions of the
blocks of one are compared with those of the chunks of the other. The
numbers each filter keeps are those of bcftools 1.24 as well, and
`tests/reference/filters/` holds every position each one keeps.

The counts of a chain are compared with pyNei's `gather_filtering_stats`,
kind by kind: pyNei gives the last filter first and popnei gives the filters
in the order of the steps, which is asserted apart, over `list(filtering)`,
because two dicts of the same pairs are equal in any order.
"""

from pathlib import Path

import numpy
import pytest
from popnei import FilteringStats, Step, open_vars, open_vcf, write_vars
from popnei.variant import Variants
from pynei import vars_from_vcf
from pynei.var_filters import (
    filter_by_maf,
    filter_by_missing_data,
    filter_by_obs_het,
    gather_filtering_stats,
)

# The variants of `many.vcf` and the size of block the tests that stop half
# way through a pass ask for: three blocks of 7 variants are 21 of them.
MANY_NUM_VARS = 500
NUM_VARS_PER_BLOCK = 7
BLOCKS_READ = 3
VARS_OF_THE_BLOCKS_READ = 21

# The method of popnei, the function of pyNei and the name of the argument
# of each of the three filters, under the kind of the filter, which is the
# name its counts have in `pass_stats.filtering`.
FILTERS = {
    "missing_data": (
        "filter_by_missing_data",
        filter_by_missing_data,
        "max_allowed_missing_rate",
    ),
    "maf": ("filter_by_maf", filter_by_maf, "max_allowed_maf"),
    "obs_het": ("filter_by_obs_het", filter_by_obs_het, "max_allowed_obs_het"),
}

# The table of "How it is verified" of `docs/specs/filters.md`: each filter
# with its threshold and how many of the 500 variants of `many.vcf` it
# keeps, which bcftools 1.24 and pyNei agree on.
THE_TABLE = [
    ("missing_data", 0.0, 26),
    ("missing_data", 0.04, 215),
    ("missing_data", 0.1, 455),
    ("maf", 0.5, 35),
    ("maf", 0.8, 384),
    ("maf", 0.95, 480),
    ("obs_het", 0.1, 22),
    ("obs_het", 0.25, 79),
    ("obs_het", 0.5, 369),
]

# The chain of "How it is verified" of the counts: the missing data filter
# at 0.04, the maf filter at 0.8 after it and the observed heterozygosity
# one at 0.5 after that, which keep 215, 163 and 106 variants of `many.vcf`.
THE_CHAIN = (("missing_data", 0.04), ("maf", 0.8), ("obs_het", 0.5))
THE_COUNTS_OF_THE_CHAIN = {
    "missing_data": FilteringStats(vars_processed=500, vars_kept=215),
    "maf": FilteringStats(vars_processed=215, vars_kept=163),
    "obs_het": FilteringStats(vars_processed=163, vars_kept=106),
}
VARS_KEPT_BY_THE_CHAIN = 106

# The thresholds that are not a number from 0 to 1, which every filter
# refuses at the call, and how each one is written in the message.
THRESHOLDS_REFUSED = [(-0.1, "-0.1"), (1.5, "1.5"), (float("nan"), "NaN")]


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which is what pyNei reads."""
    return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)


def _filtered(variants: Variants, filters) -> Variants:
    """`variants` with each of `filters`, a kind and a threshold, put on it
    in that order."""
    for kind, threshold in filters:
        getattr(variants, FILTERS[kind][0])(threshold)
    return variants


def _joined(variants: Variants):
    """The genotypes and the positions of every block of one pass, joined."""
    blocks = list(variants.iter_blocks(fields=("pos",)))
    return (
        numpy.concatenate([block.gts for block in blocks]),
        tuple(int(pos) for block in blocks for pos in block.pos),
    )


def _joined_chunks(variants):
    """The same two of the chunks of pyNei, joined."""
    chunks = list(variants.iter_vars_chunks())
    return (
        numpy.concatenate([chunk.gts.gt_values for chunk in chunks]),
        tuple(int(pos) for chunk in chunks for pos in chunk.vars_info["pos"]),
    )


@pytest.mark.parametrize(("kind", "threshold", "kept"), THE_TABLE)
def test_each_filter_keeps_the_variants_of_the_table_that_pynei_keeps(
    kind: str, threshold: float, kept: int, reference_vcf_dir: Path
) -> None:
    """One row of the table of the spec, read by both libraries.

    The genotypes and the positions have to agree variant by variant, and
    not only in how many variants there are: the two libraries work the
    number of each variant out from the same genotypes, and a filter that
    kept the wrong ones would keep as many.
    """
    ours = _filtered(_many(reference_vcf_dir), [(kind, threshold)])
    theirs = FILTERS[kind][1](vars_from_vcf(reference_vcf_dir / "many.vcf"), threshold)

    our_gts, our_positions = _joined(ours)
    their_gts, their_positions = _joined_chunks(theirs)

    assert len(our_positions) == kept
    assert our_positions == their_positions
    numpy.testing.assert_array_equal(our_gts, their_gts)

    # A second pass over the same `Variants` gives what the first gave: the
    # source is read again from its start and the steps are run again.
    again_gts, again_positions = _joined(ours)
    assert again_positions == our_positions
    numpy.testing.assert_array_equal(again_gts, our_gts)


def test_the_three_methods_return_none_and_add_their_step_in_order(
    reference_vcf_dir: Path,
) -> None:
    """What a user holds after the three calls: a `Variants` with three
    steps, each with its threshold under the name of the argument.

    Each method changes the `Variants` and returns nothing, as `list.sort`
    does, so `v2 = v1.filter_by_maf(0.95)` gives a `None` and an error at
    the next line instead of two names for one filtered object.
    """
    variants = _many(reference_vcf_dir)
    assert variants.steps == ()

    assert variants.filter_by_missing_data(0.04) is None
    assert variants.filter_by_maf(0.8) is None
    assert variants.filter_by_obs_het(0.5) is None

    assert variants.steps == (
        Step(kind="missing_data", args={"max_allowed_missing_rate": 0.04}),
        Step(kind="maf", args={"max_allowed_maf": 0.8}),
        Step(kind="obs_het", args={"max_allowed_obs_het": 0.5}),
    )


def test_the_repr_of_a_variants_names_each_filter_and_its_threshold(
    reference_vcf_dir: Path,
) -> None:
    """What a user prints in a notebook whose cells were run out of order,
    to see which filters their variants carry."""
    variants = _filtered(_many(reference_vcf_dir), THE_CHAIN)

    printed = repr(variants)

    assert "missing_data(max_allowed_missing_rate=0.04)" in printed
    assert "maf(max_allowed_maf=0.8)" in printed
    assert "obs_het(max_allowed_obs_het=0.5)" in printed
    assert str(reference_vcf_dir / "many.vcf") in printed


@pytest.mark.parametrize("kind", list(FILTERS))
@pytest.mark.parametrize(("threshold", "written"), THRESHOLDS_REFUSED)
def test_a_threshold_that_is_not_a_number_from_0_to_1_is_refused_at_the_call(
    kind: str, threshold: float, written: str, reference_vcf_dir: Path
) -> None:
    """-0.1, 1.5 and NaN, in each of the three methods.

    The number of a variant that the threshold is compared with is one
    count of the variant divided by another, so no other threshold says
    anything about which variants a user wants: pyNei takes them, and a 95
    written for 0.95 filters nothing there and says nothing. The message
    names the argument and the value, which is what shows a user that they
    wrote 95 for 0.95, and the step is not added.
    """
    method, _, argument = FILTERS[kind]
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError) as refusal:
        getattr(variants, method)(threshold)

    assert argument in str(refusal.value)
    assert written in str(refusal.value)
    assert variants.steps == ()


@pytest.mark.parametrize("kind", list(FILTERS))
def test_a_filter_with_no_threshold_is_a_type_error(
    kind: str, reference_vcf_dir: Path
) -> None:
    """Whoever filters means to filter by some rate, so no threshold is the
    natural one and none is the default, where pyNei's missing data filter
    has 0.0, which keeps the 26 variants of `many.vcf` with every genotype
    called."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError):
        getattr(variants, FILTERS[kind][0])()


def test_a_second_filter_of_one_kind_is_refused_with_the_threshold_that_is_set(
    reference_vcf_dir: Path,
) -> None:
    """Two threshold filters of one kind keep what the stricter of them
    keeps alone, so the second says that the user has lost track of what
    their `Variants` holds, which running the cell of a notebook twice
    gives.

    The message names the kind and both thresholds, the one that is set and
    the one that was written, and a filter of another kind between the two
    changes nothing.
    """
    variants = _filtered(_many(reference_vcf_dir), [("maf", 0.8)])

    with pytest.raises(ValueError) as refusal:
        variants.filter_by_maf(0.95)

    assert "maf" in str(refusal.value)
    assert "0.8" in str(refusal.value)
    assert "0.95" in str(refusal.value)
    assert variants.steps == (Step(kind="maf", args={"max_allowed_maf": 0.8}),)

    variants.filter_by_missing_data(0.04)
    with pytest.raises(ValueError, match="maf"):
        variants.filter_by_maf(0.5)


def test_the_counts_of_a_chain_are_those_of_pynei_in_the_order_of_the_steps(
    reference_vcf_dir: Path,
) -> None:
    """The three filters of "How it is verified" of the counts, at 0.04,
    0.8 and 0.5, over `many.vcf` in both libraries.

    pyNei keeps the counts in its `Variants` and gives the last filter
    first; popnei gives them with the result of the pass, in the order of
    the steps. A dict compares equal in any order, so the order is
    asserted over `list(filtering)`.
    """
    ours = _filtered(_many(reference_vcf_dir), THE_CHAIN)
    theirs = filter_by_obs_het(
        filter_by_maf(
            filter_by_missing_data(vars_from_vcf(reference_vcf_dir / "many.vcf"), 0.04),
            0.8,
        ),
        0.5,
    )

    blocks = ours.iter_blocks(fields=("pos",))
    our_positions = tuple(int(pos) for block in blocks for pos in block.pos)
    stats = blocks.pass_stats

    # pyNei counts a pass when it is made, and `gather_filtering_stats`
    # gives the counts of the last one.
    _, their_positions = _joined_chunks(theirs)
    of_pynei = gather_filtering_stats(theirs)

    assert len(our_positions) == VARS_KEPT_BY_THE_CHAIN
    assert our_positions == their_positions
    assert stats.num_vars == VARS_KEPT_BY_THE_CHAIN
    assert list(stats.filtering) == ["missing_data", "maf", "obs_het"]
    for kind, counts in THE_COUNTS_OF_THE_CHAIN.items():
        assert stats.filtering[kind] == counts
        assert of_pynei[kind].vars_processed == counts.vars_processed
        assert of_pynei[kind].vars_kept == counts.vars_kept


def test_every_pass_counts_its_own_variants_and_a_variants_with_no_filter_counts_none(
    reference_vcf_dir: Path,
) -> None:
    """A second `iter_blocks` over the same `Variants`, and one over a
    `Variants` that nothing was put on.

    Every pass builds its own filters from the steps, so the counts of the
    second pass are those of the second pass and not the double, and a
    pass with no filter has nothing to count.
    """
    variants = _many(reference_vcf_dir)

    blocks = variants.iter_blocks()
    assert sum(block.num_vars for block in blocks) == MANY_NUM_VARS
    assert blocks.pass_stats.filtering == {}

    _filtered(variants, [("missing_data", 0.04)])
    first = variants.iter_blocks()
    assert sum(block.num_vars for block in first) == 215
    second = variants.iter_blocks()
    assert sum(block.num_vars for block in second) == 215

    of_one_filter = {"missing_data": FilteringStats(vars_processed=500, vars_kept=215)}
    assert first.pass_stats.filtering == of_one_filter
    assert second.pass_stats.filtering == of_one_filter


def test_the_counts_read_while_a_pass_runs_are_of_the_blocks_it_gave(
    reference_vcf_dir: Path,
) -> None:
    """Three blocks of 7 variants of a pass with the three filters on it.

    `num_vars` is of the blocks the user got, 21 of them. What the filters
    have counted by then is of the variants the chain has read, which can be
    more, because the reader that cuts the blocks to the size the user asked
    for keeps the variants of its next block.
    """
    variants = _filtered(_many(reference_vcf_dir), THE_CHAIN)

    blocks = variants.iter_blocks(num_vars_per_block=NUM_VARS_PER_BLOCK)
    for _ in range(BLOCKS_READ):
        next(blocks)

    assert blocks.pass_stats.num_vars == VARS_OF_THE_BLOCKS_READ


def test_a_filter_added_while_a_pass_runs_holds_from_the_next_pass(
    reference_vcf_dir: Path,
) -> None:
    """A filter put on the `Variants` from inside the loop of an
    `iter_blocks`.

    The pass that runs took the steps when it started, so it gives the 500
    variants of the file and counts no filter, and the next pass has the
    filter and its counts.
    """
    variants = _many(reference_vcf_dir)
    blocks = variants.iter_blocks(num_vars_per_block=NUM_VARS_PER_BLOCK)

    given = 0
    for block in blocks:
        given += block.num_vars
        if variants.steps == ():
            variants.filter_by_missing_data(0.04)

    assert given == MANY_NUM_VARS
    assert blocks.pass_stats.num_vars == MANY_NUM_VARS
    assert blocks.pass_stats.filtering == {}

    after = variants.iter_blocks()
    assert sum(block.num_vars for block in after) == 215
    assert after.pass_stats.filtering == {
        "missing_data": FilteringStats(vars_processed=500, vars_kept=215)
    }


def test_write_vars_writes_the_variants_the_filters_kept_and_counts_them(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The test that "Its Python and TypeScript functions" of the writer in
    `docs/specs/io_vars.md` asks for: `many.vcf` with every variant given
    and the missing data filter at 0.04.

    The pass is the core's, which writes one batch for each block, so the
    count is of what was written; the file read back holds those variants
    and no step of the `Variants` that was written from.
    """
    variants = _filtered(_many(reference_vcf_dir), [("missing_data", 0.04)])
    path = tmp_path / "filtered.vars"

    written = write_vars(variants, path)

    assert written.pass_stats.num_vars == 215
    assert written.pass_stats.filtering == {
        "missing_data": FilteringStats(vars_processed=500, vars_kept=215)
    }
    read_back = open_vars(path)
    assert read_back.steps == ()
    assert sum(block.num_vars for block in read_back.iter_blocks()) == 215


def test_a_filter_over_a_vars_file_keeps_what_it_keeps_over_the_vcf(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The same filter over the two sources popnei reads.

    A filter is a reader over the reader of the source, so which variants
    it keeps depends on the genotypes and not on where they were read from,
    and the batches of the file, 100 variants each, are not the blocks the
    pass gives.
    """
    path = tmp_path / "many.vars"
    write_vars(_many(reference_vcf_dir), path, 100)

    of_the_vcf = _joined(_filtered(_many(reference_vcf_dir), [("maf", 0.5)]))
    of_the_file = _joined(_filtered(open_vars(path), [("maf", 0.5)]))

    assert of_the_file[1] == of_the_vcf[1]
    numpy.testing.assert_array_equal(of_the_file[0], of_the_vcf[0])
