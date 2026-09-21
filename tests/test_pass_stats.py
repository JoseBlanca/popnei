"""The counts of one pass, and the steps a `Variants` holds.

Every consumer of a `Variants` gives back the counts of the pass it made:
how many variants it took, and how many variants each filter of the pass was
given and kept. `docs/specs/variant.md` has the `PassStats` they come in,
`docs/specs/filters.md` the `FilteringStats` of one filter and the `Step`
that a filter is, and `docs/specs/io_vars.md` the `VarsWritten` that
`write_vars` gives.

No test here puts a filter on a `Variants`, so every `filtering` is empty
and every `steps` is `()`: what a filter counts and what it keeps are
asserted by `tests/test_filters.py`. The numbers are those of `many.vcf` of
`docs/specs/io_vcf.md`, 500 variants of 50 individuals read with every
variant given.
"""

import dataclasses
from pathlib import Path

import pytest
from popnei import (
    FilteringStats,
    PassStats,
    Step,
    VarsWritten,
    open_vars,
    open_vcf,
    write_vars,
)
from popnei.variant import Variants, _pass_stats_of

# The VCFs that `tests/reference/vcf/make_reference.py` writes, which
# `conftest.py` gives the tests as `reference_vcf_dir` and which the module
# scoped fixture below reads by this path.
REFERENCE_VCF_DIR = Path(__file__).parent / "reference" / "vcf"

# The variants of `many.vcf`, read with every variant given, and the size of
# block the tests that stop half way through a pass use: three blocks of 7
# variants are 21 of the 500.
MANY_NUM_VARS = 500
NUM_VARS_PER_BLOCK = 7
BLOCKS_READ = 3
VARS_OF_THE_BLOCKS_READ = 21

# How many variants a batch of the vars file the tests write holds, which is
# not the size of the blocks any of them asks for: the counts are of the
# blocks a user got and not of what the reader took from the file.
VARS_NUM_VARS_PER_BLOCK = 100


@pytest.fixture(scope="module")
def many_vars(tmp_path_factory) -> Path:
    """`many.vcf`, every variant of it, as a vars file of batches of 100.

    It is written once for the module: the passes over it read it and do
    not change it.
    """
    path = tmp_path_factory.mktemp("vars") / "many.vars"
    variants = open_vcf(REFERENCE_VCF_DIR / "many.vcf", only_passed=False)
    write_vars(variants, path, VARS_NUM_VARS_PER_BLOCK)
    return path


@pytest.fixture(params=["a vcf", "a vars file"])
def many(request, reference_vcf_dir: Path, many_vars: Path) -> Variants:
    """The 500 variants of `many.vcf`, from the VCF and from a vars file.

    Every source gives the same counts: what a pass counts is the variants
    of the blocks it gave, whatever file they were read from.
    """
    if request.param == "a vcf":
        return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)
    return open_vars(many_vars)


def test_a_whole_iter_blocks_counts_every_variant_it_gave(many: Variants) -> None:
    """The 500 variants of the source, and no filter to count."""
    blocks = many.iter_blocks()

    assert sum(block.num_vars for block in blocks) == MANY_NUM_VARS
    assert blocks.pass_stats == PassStats(num_vars=MANY_NUM_VARS, filtering={})


def test_the_counts_of_a_pass_that_is_not_over_are_of_the_blocks_it_gave(
    many: Variants,
) -> None:
    """Three blocks of 7 variants, read out of a source of 500.

    The counts are read while the pass runs, which is what a user does in
    the loop of an `iter_blocks`, and they hold the 21 variants of the
    three blocks and not the variants the reader has taken out of the file:
    a batch of the vars file is 100 variants.
    """
    blocks = many.iter_blocks(num_vars_per_block=NUM_VARS_PER_BLOCK)

    assert blocks.pass_stats == PassStats(num_vars=0, filtering={})
    for _ in range(BLOCKS_READ):
        next(blocks)
    assert blocks.pass_stats == PassStats(
        num_vars=VARS_OF_THE_BLOCKS_READ, filtering={}
    )


def test_every_pass_counts_its_own_variants(many: Variants) -> None:
    """A second `iter_blocks` gives 500 again and not 1000.

    Every call starts a pass of its own, with counts of its own, so nothing
    of the pass before it is added to them.
    """
    first = many.iter_blocks()
    assert sum(block.num_vars for block in first) == MANY_NUM_VARS

    second = many.iter_blocks()
    assert sum(block.num_vars for block in second) == MANY_NUM_VARS
    assert second.pass_stats == PassStats(num_vars=MANY_NUM_VARS, filtering={})
    assert first.pass_stats == PassStats(num_vars=MANY_NUM_VARS, filtering={})


def test_write_vars_gives_the_counts_of_the_pass_it_made(
    many: Variants, tmp_path: Path
) -> None:
    """The whole source read into a vars file, and the 500 variants of it.

    The pass is the core's, which writes one batch for each block, so the
    count is what was written and not what a loop of Python saw.
    """
    path = tmp_path / "written.vars"

    written = write_vars(many, path, VARS_NUM_VARS_PER_BLOCK)

    assert isinstance(written, VarsWritten)
    assert written.pass_stats == PassStats(num_vars=MANY_NUM_VARS, filtering={})
    assert sum(block.num_vars for block in open_vars(path).iter_blocks()) == (
        MANY_NUM_VARS
    )


def test_a_source_with_no_variants_counts_none(write_vcf, tmp_path: Path) -> None:
    """A VCF whose header names three individuals and that has no variant.

    It is not an error, so the pass is a pass like any other and its count
    is 0.
    """
    variants = open_vcf(write_vcf([]))

    blocks = variants.iter_blocks()
    assert list(blocks) == []
    assert blocks.pass_stats == PassStats(num_vars=0, filtering={})

    written = write_vars(variants, tmp_path / "empty.vars")
    assert written.pass_stats == PassStats(num_vars=0, filtering={})


def test_the_block_a_pass_lost_with_an_error_is_not_among_its_variants(
    write_vcf,
) -> None:
    """Three variants read in blocks of one, and a fourth line popnei
    refuses because it gives one individual four alleles.

    The block the error happened in never reached the user, so the count is
    of the three blocks they got and not of the four variants the file
    holds.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0/1/1\t1/1",
        ]
    )
    blocks = open_vcf(path).iter_blocks(num_vars_per_block=1)

    given = []
    with pytest.raises(ValueError, match="ind2"):
        for block in blocks:
            given.append(block.num_vars)

    assert given == [1, 1, 1]
    assert blocks.pass_stats == PassStats(num_vars=3, filtering={})


def test_the_counts_of_a_pass_are_frozen_dataclasses() -> None:
    """`PassStats`, `FilteringStats` and `VarsWritten`: a result of popnei
    is not changed.

    A user who reads the counts of a pass into a report of their own gets
    numbers that nothing can write over, and two results of the same pass
    are equal.
    """
    stats = PassStats(num_vars=3, filtering={"maf": FilteringStats(5, 3)})
    written = VarsWritten(pass_stats=stats)

    assert dataclasses.is_dataclass(PassStats)
    assert dataclasses.is_dataclass(FilteringStats)
    assert dataclasses.is_dataclass(VarsWritten)
    assert stats == PassStats(num_vars=3, filtering={"maf": FilteringStats(5, 3)})
    assert FilteringStats(vars_processed=5, vars_kept=3) == FilteringStats(5, 3)
    assert written == VarsWritten(pass_stats=stats)
    with pytest.raises(dataclasses.FrozenInstanceError):
        stats.num_vars = 4
    with pytest.raises(dataclasses.FrozenInstanceError):
        FilteringStats(5, 3).vars_kept = 5
    with pytest.raises(dataclasses.FrozenInstanceError):
        written.pass_stats = stats


def test_the_counts_of_the_filters_come_in_the_order_of_the_steps() -> None:
    """The chain of readers gives the outermost filter first, and a user
    reads the filters in the order in which they were put on the
    `Variants`.

    The numbers are those of `docs/specs/filters.md`, the missing data
    filter at 0.04 and the maf filter at 0.8 on `many.vcf`: 500 variants
    given and 215 kept, and then 215 given and 163 kept. The chain has the
    maf filter first, because it is the outermost, and `filtering` has the
    missing data one first.
    """
    stats = _pass_stats_of((163, [("maf", 215, 163), ("missing_data", 500, 215)]))

    assert list(stats.filtering) == ["missing_data", "maf"]
    assert stats.filtering["missing_data"] == FilteringStats(
        vars_processed=500, vars_kept=215
    )
    assert stats.filtering["maf"] == FilteringStats(vars_processed=215, vars_kept=163)
    assert stats.num_vars == 163


def test_iter_blocks_refuses_its_arguments_at_the_call_and_not_at_a_block(
    reference_vcf_dir: Path,
) -> None:
    """A field that is of nothing and a block of no variants.

    The call itself raises, with no block asked for: a user who wrote a
    name wrong is told where they wrote it, and not at the `for` that comes
    later.
    """
    variants = open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)

    with pytest.raises(ValueError, match="depth"):
        variants.iter_blocks(fields=("chrom", "depth"))
    with pytest.raises(ValueError, match="0 variants"):
        variants.iter_blocks(num_vars_per_block=0)
    with pytest.raises(TypeError, match="sequence"):
        variants.iter_blocks(fields="alleles")


def test_a_variants_just_opened_has_no_step(many: Variants) -> None:
    """The steps of a source that nothing was asked of: none of them.

    `steps` is a tuple, so what a user reads is not a list of popnei's that
    they could add to: a step is put on a `Variants` by a method of it.
    """
    assert many.steps == ()
    assert isinstance(many.steps, tuple)


def test_a_step_is_a_frozen_dataclass_of_a_kind_and_its_arguments() -> None:
    """`Step`, which `variants.steps` holds one of for each step.

    `args` is a dict of the name of each argument, as the user writes it,
    to its value, so that the steps of the later filters, which take other
    arguments, fit in it.
    """
    step = Step(kind="maf", args={"max_allowed_maf": 0.95})

    assert dataclasses.is_dataclass(Step)
    assert step == Step(kind="maf", args={"max_allowed_maf": 0.95})
    assert step.kind == "maf"
    assert step.args == {"max_allowed_maf": 0.95}
    with pytest.raises(dataclasses.FrozenInstanceError):
        step.kind = "obs_het"


def test_the_repr_of_a_variants_names_its_source_and_says_it_has_no_step(
    reference_vcf_dir: Path, many_vars: Path
) -> None:
    """What a user prints in a notebook to see what a handle holds.

    A notebook whose cells were run out of order is where a user needs to
    read which steps are on a `Variants`, so the `repr` says that there are
    none as plainly as it would name them.
    """
    vcf_path = reference_vcf_dir / "many.vcf"

    printed = repr(open_vcf(vcf_path, only_passed=False))
    assert str(vcf_path) in printed
    assert "no steps" in printed

    printed = repr(open_vars(many_vars))
    assert str(many_vars) in printed
    assert "no steps" in printed
