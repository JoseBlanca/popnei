"""The filter of individuals from Python: whose genotypes come out of a
pass, in which order, and what a `Variants` carries once it is put on.

`docs/specs/filters.md` has the filter under "The filter of individuals".
The comparison it asks for is with pyNei at commit ef0ca6e, which
`pyproject.toml` names: `many.vcf` of `docs/specs/io_vcf.md`, 500 variants
of 50 diploid individuals, is read by both libraries, popnei with
`only_passed=False` because pyNei gives every variant whatever its FILTER
says, and `ind05`, `ind00` and `ind49` are kept in both.

The genotypes are compared column by name and not column by column: popnei
keeps the individuals in the order the user named them and pyNei's
`filter_samples` keeps them in the order of the source, so the two arrays
hold the same three columns in a different order. The numbers a filter of
individuals with the missing data filter after it keeps are those of
bcftools 1.24 as well, which "How it is verified" of the filter gives.
"""

import copy
from pathlib import Path

import numpy
import pytest
from popnei import FilteringStats, Step, open_vcf
from popnei.variant import Variants
from pynei import vars_from_vcf
from pynei.var_filters import filter_by_missing_data, filter_samples

# The three individuals of "How it is verified" of the filter, in the order
# a user names them, which is not the order `many.vcf` has them in.
THE_THREE = ("ind05", "ind00", "ind49")
THE_THREE_OF_THE_SOURCE = ("ind00", "ind05", "ind49")

# The 500 variants of `many.vcf`, and what the missing data filter at 0
# keeps after the three individuals are taken: 423 variants, at these five
# positions first. The same filter over the 50 individuals keeps 26.
MANY_NUM_VARS = 500
KEPT_OF_THE_THREE = 423
FIRST_POSITIONS_OF_THE_THREE = (1000, 1037, 1074, 1111, 1148)
KEPT_OF_THE_FIFTY = 26

# The genotypes of `ind05`, `ind00` and `ind49`, in that order, at the first
# variant the missing data filter keeps and at the third, which are `1|1`,
# `1/1`, `1/1` and `0/1`, `2|1`, `1|2` in the VCF. The third variant is the
# one that tells the order of the argument from the order of the source.
GTS_AT_1000 = [[1, 1], [1, 1], [1, 1]]
GTS_AT_1074 = [[0, 1], [2, 1], [1, 2]]


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which is what pyNei reads."""
    return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)


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


def test_the_genotypes_kept_are_those_pynei_keeps_for_the_same_individuals(
    reference_vcf_dir: Path,
) -> None:
    """`many.vcf` with `ind05`, `ind00` and `ind49` kept, in both libraries.

    Each column of popnei's genotypes is compared with the column of the
    same individual in pyNei's, which is at another place: popnei gives the
    three in the order of the argument and pyNei in the order of the source.
    Every variant of the file comes out, since the filter takes no variant
    away, and the positions have to agree variant by variant.
    """
    ours = _many(reference_vcf_dir)
    ours.filter_individuals(THE_THREE)
    theirs = filter_samples(
        vars_from_vcf(reference_vcf_dir / "many.vcf"), samples=list(THE_THREE)
    )

    our_gts, our_positions = _joined(ours)
    their_gts, their_positions = _joined_chunks(theirs)
    their_names = list(theirs.samples)

    assert ours.individuals == THE_THREE
    assert tuple(their_names) == THE_THREE_OF_THE_SOURCE
    assert len(our_positions) == MANY_NUM_VARS
    assert our_positions == their_positions
    assert our_gts.shape == (MANY_NUM_VARS, 3, 2)
    for column, name in enumerate(ours.individuals):
        numpy.testing.assert_array_equal(
            our_gts[:, column, :], their_gts[:, their_names.index(name), :]
        )


def test_the_missing_data_filter_after_it_keeps_the_423_variants_and_counts_them(
    reference_vcf_dir: Path,
) -> None:
    """The three individuals with the missing data filter at 0 after them.

    Of the 500 variants, 423 have every genotype of the three called, where
    only 26 have every genotype of the 50 called: the filter after the step
    divides by the kept individuals alone. The filter of individuals takes
    no variant away and has no counts, so the counts of the pass hold the
    missing data filter alone, with the 500 it was given and the 423 it
    kept.
    """
    ours = _many(reference_vcf_dir)
    ours.filter_individuals(THE_THREE)
    ours.filter_by_missing_data(0)
    theirs = filter_by_missing_data(
        filter_samples(
            vars_from_vcf(reference_vcf_dir / "many.vcf"), samples=list(THE_THREE)
        ),
        0,
    )

    blocks = ours.iter_blocks(fields=("pos",))
    kept = list(blocks)
    our_positions = tuple(int(pos) for block in kept for pos in block.pos)
    stats = blocks.pass_stats
    _, their_positions = _joined_chunks(theirs)

    assert len(our_positions) == KEPT_OF_THE_THREE
    assert our_positions[:5] == FIRST_POSITIONS_OF_THE_THREE
    assert our_positions == their_positions
    assert stats.num_vars == KEPT_OF_THE_THREE
    assert stats.filtering == {
        "missing_data": FilteringStats(
            vars_processed=MANY_NUM_VARS, vars_kept=KEPT_OF_THE_THREE
        )
    }


def test_the_genotypes_come_in_the_order_of_the_argument(
    reference_vcf_dir: Path,
) -> None:
    """The two variants of the spec, at the positions 1000 and 1074, of
    `ind05`, `ind00` and `ind49`.

    At 1074 the three genotypes differ, so the row says which column holds
    which individual: keeping them in the order of the source, as pyNei
    does, would put `ind00` first.
    """
    variants = _many(reference_vcf_dir)
    variants.filter_individuals(THE_THREE)
    variants.filter_by_missing_data(0)

    gts, positions = _joined(variants)

    assert positions[0] == 1000
    numpy.testing.assert_array_equal(gts[0], GTS_AT_1000)
    numpy.testing.assert_array_equal(gts[positions.index(1074)], GTS_AT_1074)


def test_a_filter_before_it_counts_over_every_individual_of_the_source(
    reference_vcf_dir: Path,
) -> None:
    """The missing data filter at 0 put on before the three individuals.

    A step sees what the steps before it gave, so this filter divides by the
    50 individuals of the source and keeps the 26 variants with every
    genotype called, where the same filter after the step keeps 423.
    """
    variants = _many(reference_vcf_dir)
    variants.filter_by_missing_data(0)
    variants.filter_individuals(THE_THREE)

    blocks = variants.iter_blocks(fields=("pos",))
    kept = list(blocks)

    assert sum(block.num_vars for block in kept) == KEPT_OF_THE_FIFTY
    assert blocks.pass_stats.filtering == {
        "missing_data": FilteringStats(
            vars_processed=MANY_NUM_VARS, vars_kept=KEPT_OF_THE_FIFTY
        )
    }
    assert kept[0].gts.shape[1] == 3


def test_the_method_returns_none_and_adds_its_step_with_the_names_as_a_tuple(
    reference_vcf_dir: Path,
) -> None:
    """What a user holds after the call: a `Variants` with one step, whose
    arguments are the names they wrote, in their order.

    The method changes the `Variants` and returns nothing, as the three
    threshold filters do, so `v2 = v1.filter_individuals(names)` gives a
    `None` and an error at the next line instead of two names for one
    filtered object. The filter of individuals takes no variant away, so it
    has no entry in the counts of a pass.
    """
    variants = _many(reference_vcf_dir)
    assert variants.steps == ()

    assert variants.filter_individuals(list(THE_THREE)) is None

    assert variants.steps == (
        Step(kind="individuals", args={"individuals": THE_THREE}),
    )
    assert variants.iter_blocks().pass_stats.filtering == {}


def test_the_individuals_are_the_kept_ones_before_a_pass_and_after_one(
    reference_vcf_dir: Path,
) -> None:
    """`individuals` and `num_individuals` read before the first pass and
    after a whole one.

    They are what the next pass gives and not what the header holds, and
    nothing of a pass changes them: the names are a tuple, as they are for a
    `Variants` that nothing was put on.
    """
    variants = _many(reference_vcf_dir)
    assert variants.num_individuals == 50
    assert variants.individuals[:2] == ("ind00", "ind01")

    variants.filter_individuals(THE_THREE)

    assert variants.individuals == THE_THREE
    assert variants.num_individuals == 3
    assert isinstance(variants.individuals, tuple)

    for block in variants.iter_blocks():
        assert block.gts.shape[1] == 3

    assert variants.individuals == THE_THREE
    assert variants.num_individuals == 3


def test_two_handles_over_one_file_each_keep_their_own_individuals(
    reference_vcf_dir: Path,
) -> None:
    """The two populations of pyNei's
    `test_several_vars_can_share_one_filtered_source`, which puts two
    filters of individuals over one `Variants`.

    A second filter of individuals on one `Variants` is refused here, so a
    user who wants two sets of individuals opens the source twice, which
    reads the header and nothing else. Each handle keeps its own.
    """
    first = _many(reference_vcf_dir)
    second = _many(reference_vcf_dir)
    first.filter_individuals(("ind00", "ind01"))
    second.filter_individuals(("ind02", "ind03", "ind04"))

    of_the_first, _ = _joined(first)
    of_the_second, _ = _joined(second)
    whole, _ = _joined(_many(reference_vcf_dir))

    assert first.individuals == ("ind00", "ind01")
    assert second.individuals == ("ind02", "ind03", "ind04")
    numpy.testing.assert_array_equal(of_the_first, whole[:, 0:2, :])
    numpy.testing.assert_array_equal(of_the_second, whole[:, 2:5, :])


def test_the_repr_of_a_variants_names_the_step_and_the_individuals_kept(
    reference_vcf_dir: Path,
) -> None:
    """What a user prints in a notebook whose cells were run out of order,
    to see which individuals their variants carry."""
    variants = _many(reference_vcf_dir)
    variants.filter_individuals(THE_THREE)

    printed = repr(variants)

    assert f"individuals(individuals={THE_THREE!r})" in printed
    assert str(reference_vcf_dir / "many.vcf") in printed


def test_a_name_that_is_not_an_individual_of_the_source_is_refused_at_the_call(
    reference_vcf_dir: Path,
) -> None:
    """A name the source does not have, beside one it has.

    pyNei drops it in silence and gives a `Variants` of one individual, so
    a typed name is a result over the wrong individuals there. The message
    names what the user wrote, the step is not added and the individuals are
    those of the source.
    """
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError) as refusal:
        variants.filter_individuals(("ind05", "nope"))

    assert "nope" in str(refusal.value)
    assert variants.steps == ()
    assert variants.num_individuals == 50


def test_a_name_that_is_there_twice_is_refused_at_the_call(
    reference_vcf_dir: Path,
) -> None:
    """One individual named twice.

    Two columns of the genotypes of one individual are one individual for
    everything that reads them, and every count over them would hold it
    twice. pyNei keeps the individual once. The message names the name.
    """
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError) as refusal:
        variants.filter_individuals(("ind05", "ind00", "ind05"))

    assert "ind05" in str(refusal.value)
    assert variants.steps == ()


def test_a_filter_of_no_individual_is_refused_at_the_call(
    reference_vcf_dir: Path,
) -> None:
    """No name at all, which would leave variants of nobody: every source of
    popnei holds one individual at least."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError):
        variants.filter_individuals(())

    assert variants.steps == ()


def test_a_second_filter_of_individuals_is_refused_with_the_kind(
    reference_vcf_dir: Path,
) -> None:
    """Two filters of individuals on one `Variants`.

    Two lists keep the individuals that are in both, which is one list, so
    the second says that the user has lost track of the individuals their
    variants carry, which running the cell of a notebook twice gives. pyNei
    takes it. After the refusal the steps are as they were, and a threshold
    filter between the two changes nothing.
    """
    variants = _many(reference_vcf_dir)
    variants.filter_individuals(THE_THREE)

    with pytest.raises(ValueError) as refusal:
        variants.filter_individuals(("ind01", "ind02"))

    assert "individuals" in str(refusal.value)
    assert variants.steps == (
        Step(kind="individuals", args={"individuals": THE_THREE}),
    )
    assert variants.individuals == THE_THREE

    variants.filter_by_maf(0.95)
    with pytest.raises(ValueError, match="individuals"):
        variants.filter_individuals(("ind01",))


def test_a_shallow_copy_shares_the_steps_and_the_kept_individuals(
    reference_vcf_dir: Path,
) -> None:
    """`copy.copy(variants)` after a filter of individuals.

    A shallow copy is a second handle over the same source and the same
    steps, and the individuals a handle gives are the ones its steps keep,
    so a filter of individuals put on either of them is on both and both
    give the kept names, in the order they were named.
    """
    variants = _many(reference_vcf_dir)
    twin = copy.copy(variants)

    variants.filter_individuals(THE_THREE)

    assert twin.steps == variants.steps
    assert twin.individuals == THE_THREE
    assert twin.num_individuals == len(THE_THREE)


def test_one_name_written_as_a_string_is_a_type_error(
    reference_vcf_dir: Path,
) -> None:
    """`filter_individuals("ind05")` instead of `("ind05",)`.

    A string is a sequence of its letters, so the call would ask for the
    individuals `i`, `n`, `d`, `0` and `5`, and the user would read that `i`
    is not an individual of the variants. The message says what to write
    instead, as `iter_blocks` does for one field.
    """
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError) as refusal:
        variants.filter_individuals("ind05")

    assert "ind05" in str(refusal.value)
    assert variants.steps == ()


def test_what_is_no_sequence_of_names_is_a_type_error_that_names_the_argument(
    reference_vcf_dir: Path,
) -> None:
    """`filter_individuals(5)` and `filter_individuals(["ind00", 3])`.

    Python says `'int' object is not iterable` of the first and pyo3 says
    `'int' object is not an instance of 'str'` of the second, and neither
    message names the argument or the call the user wrote. The package
    checks both, as the TypeScript one does.
    """
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError) as of_no_sequence:
        variants.filter_individuals(5)  # type: ignore[arg-type]

    assert "individuals" in str(of_no_sequence.value)
    assert "5" in str(of_no_sequence.value)
    assert "int" in str(of_no_sequence.value)

    with pytest.raises(TypeError) as of_a_name_that_is_no_name:
        variants.filter_individuals(["ind00", 3])  # type: ignore[list-item]

    assert "individuals" in str(of_a_name_that_is_no_name.value)
    assert "3" in str(of_a_name_that_is_no_name.value)
    assert "int" in str(of_a_name_that_is_no_name.value)
    assert variants.steps == ()
