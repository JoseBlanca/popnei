"""The filter of the first n variants from Python: which variants it keeps,
the counts of a pass it ended, and what it refuses at the call.

`docs/specs/filters.md` has the filter under "The filter that keeps the
first n variants". pyNei has no such filter, so the variants are those of
bcftools 1.24 on `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
diploid individuals read with every variant given, as "How it is verified"
of the spec gives them: `bcftools view -H many.vcf | head -n 10` and, with
the MAF filter of 0.8 before the first 10, `bcftools view -H -Q 0.8:major
many.vcf | head -n 10`, run on 5 October 2026.
"""

from pathlib import Path

import pandas
import pytest
from popnei import (
    FilteringStats,
    PassStats,
    Step,
    calc_gwas,
    calc_kinship,
    calc_ld_and_dist_per_pop,
    calc_pairwise_kosman_dists,
    calc_per_individual_stats,
    calc_per_var_distribs,
    calc_pop_dists,
    calc_pop_diversity,
    calc_rogers_huff_r2_matrix,
    calc_var_density,
    do_pca_from_variants,
    do_pcoa_from_variants,
    open_vcf,
    write_vars,
    write_vcf,
)
from popnei.variant import Variants

# The first ten variants of `many.vcf`, all of chr1, by bcftools.
FIRST_TEN = [1000, 1037, 1074, 1111, 1148, 1185, 1222, 1259, 1296, 1333]

# The first ten of those the MAF filter of 0.8 keeps, by bcftools.
FIRST_TEN_AFTER_THE_MAF_FILTER = [
    1037,
    1074,
    1111,
    1148,
    1222,
    1259,
    1296,
    1333,
    1370,
    1407,
]

MANY_NUM_VARS = 500


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which is what bcftools reads."""
    return open_vcf(reference_vcf_dir / "many.vcf", only_passed=False)


def _kept(variants: Variants, num_vars_per_block: int | None = None):
    """The chromosome and the position of every variant of one pass, and
    the counts of that pass."""
    blocks = variants.iter_blocks(
        fields=("chrom", "pos"), num_vars_per_block=num_vars_per_block
    )
    kept = [
        (chrom, int(pos))
        for block in blocks
        for chrom, pos in zip(block.chrom, block.pos, strict=True)
    ]
    return kept, blocks.pass_stats


@pytest.mark.parametrize("num_vars_per_block", [7, None])
def test_filter_first_n_keeps_the_first_ten_of_bcftools_and_says_it_stopped(
    reference_vcf_dir: Path, num_vars_per_block: int | None
) -> None:
    """The ten positions, whatever the size of the blocks, and counts of
    ten variants that say the filter ended the pass."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    kept, pass_stats = _kept(variants, num_vars_per_block)

    assert kept == [("chr1", pos) for pos in FIRST_TEN]
    assert pass_stats.num_vars == 10
    assert pass_stats.stopped_early is True
    assert pass_stats.filtering["first_n"].vars_kept == 10


def test_filter_first_n_counts_the_blocks_it_took_in_blocks_of_7(
    reference_vcf_dir: Path,
) -> None:
    """With blocks of 7 the filter is given two blocks, 14 variants, and
    keeps 10; the MAF filter before it counts the same two blocks, of which
    it keeps 12, which the spec has from bcftools over the first 14 lines."""
    alone = _many(reference_vcf_dir)
    alone.filter_first_n(10)
    after_the_maf_filter = _many(reference_vcf_dir)
    after_the_maf_filter.filter_by_maf(0.8)
    after_the_maf_filter.filter_first_n(10)

    _, of_it_alone = _kept(alone, 7)
    kept, of_the_two = _kept(after_the_maf_filter, 7)

    assert of_it_alone == PassStats(
        num_vars=10,
        filtering={"first_n": FilteringStats(vars_processed=14, vars_kept=10)},
        stopped_early=True,
    )
    assert kept == [("chr1", pos) for pos in FIRST_TEN_AFTER_THE_MAF_FILTER]
    assert of_the_two == PassStats(
        num_vars=10,
        filtering={
            "maf": FilteringStats(vars_processed=14, vars_kept=12),
            "first_n": FilteringStats(vars_processed=12, vars_kept=10),
        },
        stopped_early=True,
    )


def test_filter_first_n_of_more_than_the_file_holds_keeps_them_all_and_did_not_stop(
    reference_vcf_dir: Path,
) -> None:
    """A `num_vars` of 1000 on the 500 variants: every one is kept, the
    pass reads the whole source, and the counts say it was not ended."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(1000)

    kept, pass_stats = _kept(variants, 7)

    assert len(kept) == MANY_NUM_VARS
    assert pass_stats == PassStats(
        num_vars=MANY_NUM_VARS,
        filtering={
            "first_n": FilteringStats(
                vars_processed=MANY_NUM_VARS, vars_kept=MANY_NUM_VARS
            )
        },
        stopped_early=False,
    )


def test_filter_first_n_puts_stopped_early_in_the_repr_of_the_counts(
    reference_vcf_dir: Path,
) -> None:
    """A user who prints the counts sees that the pass was ended."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    _, pass_stats = _kept(variants, 7)

    assert "stopped_early=True" in repr(pass_stats)


def test_filter_first_n_is_a_step_of_its_kind_and_num_vars(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(1000)

    assert variants.steps == (Step(kind="first_n", args={"num_vars": 1000}),)
    assert "first_n(num_vars=1000)" in repr(variants)


def test_filter_first_n_gives_a_pca_whose_second_pass_reads_the_same_variants(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """`do_pca_from_variants` reads the source twice when it gives the
    weights of the variants. `many.vcf` holds variants of three alleles, so
    the dosages count every allele that is not the major one. Over the first 50 it gives what it gives over
    a file of those 50 alone, which `write_vcf` writes through the same
    filter: a second pass that read other variants would give other
    weights."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(50)
    the_first_50 = tmp_path / "first_50.vcf"
    written = write_vcf(variants, the_first_50)

    over_the_filter = do_pca_from_variants(variants, transform_to_biallelic=True)
    over_the_file = do_pca_from_variants(
        open_vcf(the_first_50, only_passed=False), transform_to_biallelic=True
    )

    assert written.pass_stats.num_vars == 50
    assert over_the_filter.pass_stats == PassStats(
        num_vars=50,
        filtering={
            "first_n": FilteringStats(vars_processed=MANY_NUM_VARS, vars_kept=50)
        },
        stopped_early=True,
    )
    pandas.testing.assert_frame_equal(
        over_the_filter.projections, over_the_file.projections
    )
    pandas.testing.assert_frame_equal(
        over_the_filter.princomps, over_the_file.princomps
    )


@pytest.mark.parametrize(
    ("num_vars", "in_the_message"),
    [(0, "`num_vars` is 0, .* 1 or more"), (-1, "`num_vars` is -1")],
)
def test_filter_first_n_refuses_a_num_vars_below_1(
    reference_vcf_dir: Path, num_vars: int, in_the_message: str
) -> None:
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError, match=in_the_message):
        variants.filter_first_n(num_vars)
    assert variants.steps == ()


def test_filter_first_n_refuses_a_num_vars_that_64_bits_do_not_hold(
    reference_vcf_dir: Path,
) -> None:
    """The core holds the n in 64 bits on every platform, so the limit a user
    is told is 2^64 - 1, and not what the machine counts, which in pyodide
    is 2^32 - 1."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(
        ValueError, match=f"`num_vars` is {2**64}, .* 18446744073709551615"
    ):
        variants.filter_first_n(2**64)
    assert variants.steps == ()


def test_filter_first_n_takes_the_largest_num_vars_of_64_bits(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(2**64 - 1)

    assert variants.steps == (Step(kind="first_n", args={"num_vars": 2**64 - 1}),)


@pytest.mark.parametrize("num_vars", [1.5, True, "10", None])
def test_filter_first_n_refuses_what_is_no_whole_number_with_a_type_error(
    reference_vcf_dir: Path, num_vars: object
) -> None:
    """A float and a truth value are a `TypeError` that names the argument,
    as for every argument of popnei that takes a whole number."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError, match="`num_vars`"):
        variants.filter_first_n(num_vars)
    assert variants.steps == ()


def test_filter_first_n_refuses_a_second_filter_of_its_kind(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    with pytest.raises(
        ValueError,
        match="first_n already, of the first 10 variants, and a second filter of "
        "that kind, of the first 20,",
    ):
        variants.filter_first_n(20)
    assert variants.steps == (Step(kind="first_n", args={"num_vars": 10}),)


# Each step that takes variants out, as a call on a `Variants`, and the kind
# the refusal names.
_STEPS_THAT_TAKE_VARIANTS_OUT = {
    "missing_data": lambda variants: variants.filter_by_missing_data(0.1),
    "maf": lambda variants: variants.filter_by_maf(0.8),
    "obs_het": lambda variants: variants.filter_by_obs_het(0.5),
    "ld": lambda variants: variants.filter_by_ld(0.1, 10000),
    "regions": lambda variants: variants.filter_by_regions(_REGIONS_BED),
    "excluded_regions": lambda variants: variants.filter_by_regions(
        _REGIONS_BED, exclude=True
    ),
}
_REGIONS_BED = Path(__file__).parent / "reference" / "filters" / "regions.bed"


@pytest.mark.parametrize("kind", list(_STEPS_THAT_TAKE_VARIANTS_OUT))
def test_filter_first_n_refuses_a_step_that_takes_variants_out_after_it(
    reference_vcf_dir: Path, kind: str
) -> None:
    """The step is refused at the call, with its kind and the filter of the
    first n named, and the steps are as they were."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    with pytest.raises(ValueError) as refusal:
        _STEPS_THAT_TAKE_VARIANTS_OUT[kind](variants)

    assert f"a filter by {kind} after it" in str(refusal.value)
    assert "first_n already" in str(refusal.value)
    assert [step.kind for step in variants.steps] == ["first_n"]


def test_filter_first_n_refuses_a_second_maf_filter_after_it_as_a_second_of_its_kind(
    reference_vcf_dir: Path,
) -> None:
    """A MAF filter after a MAF filter and the filter of the first n breaks
    both rules, and is refused as a second filter of its kind, with both
    thresholds, as in TypeScript."""
    variants = _many(reference_vcf_dir)
    variants.filter_by_maf(0.9)
    variants.filter_first_n(10)

    with pytest.raises(ValueError) as refusal:
        variants.filter_by_maf(0.8)

    assert "filtered by maf already, with a threshold of 0.9" in str(refusal.value)
    assert "whose threshold is 0.8" in str(refusal.value)
    assert [step.kind for step in variants.steps] == ["maf", "first_n"]


def test_filter_first_n_accepts_the_filter_of_individuals_after_it(
    reference_vcf_dir: Path,
) -> None:
    """The filter of individuals takes out no variant, so the first ten are
    still ten, of the two individuals kept."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)
    variants.filter_individuals(("ind05", "ind00"))

    blocks = variants.iter_blocks(fields=("chrom", "pos"))
    gts_shapes = [block.gts.shape for block in blocks]

    assert [step.kind for step in variants.steps] == ["first_n", "individuals"]
    assert gts_shapes == [(10, 2, 2)]
    assert blocks.pass_stats.stopped_early is True


def _every_consumer(variants: Variants, tmp_path: Path) -> dict[str, PassStats]:
    """The counts each consumer of the package gives for a pass over
    `variants`."""
    individuals = variants.individuals
    pops = {"pop1": individuals[:25], "pop2": individuals[25:]}
    trait = pandas.Series(
        [float(index % 7) for index in range(len(individuals))],
        index=list(individuals),
    )
    blocks = variants.iter_blocks()
    for _ in blocks:
        pass
    kosman = calc_pairwise_kosman_dists(variants)
    return {
        "iter_blocks": blocks.pass_stats,
        "write_vars": write_vars(variants, tmp_path / "written.vars").pass_stats,
        "write_vcf": write_vcf(variants, tmp_path / "written.vcf").pass_stats,
        "calc_gwas": calc_gwas(
            variants, trait, "continuous", transform_to_biallelic=True
        ).pass_stats,
        "calc_kinship": calc_kinship(variants, transform_to_biallelic=True).pass_stats,
        "calc_ld_and_dist_per_pop": calc_ld_and_dist_per_pop(variants).pass_stats,
        "calc_pairwise_kosman_dists": kosman.pass_stats,
        "calc_per_individual_stats": calc_per_individual_stats(variants).pass_stats,
        "calc_per_var_distribs": calc_per_var_distribs(variants).pass_stats,
        "calc_pop_dists": calc_pop_dists(variants, pops, None).pass_stats,
        "calc_pop_diversity": calc_pop_diversity(variants, pops).pass_stats,
        "calc_rogers_huff_r2_matrix": calc_rogers_huff_r2_matrix(variants).pass_stats,
        "calc_var_density": calc_var_density(variants, 1000).pass_stats,
        "do_pca_from_variants": do_pca_from_variants(
            variants, transform_to_biallelic=True
        ).pass_stats,
        "do_pcoa_from_variants": do_pcoa_from_variants(
            variants, correct_by_lingoes=True
        ).pass_stats,
    }


def test_filter_first_n_is_said_to_have_stopped_by_every_consumer(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Each consumer builds its counts through the one function of the
    binding crate that reads `stopped_early`: one that did not would give
    false here."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(100)

    of_each = _every_consumer(variants, tmp_path)

    assert {name: counts.stopped_early for name, counts in of_each.items()} == (
        dict.fromkeys(of_each, True)
    )
    assert {
        name: (counts.num_vars, counts.filtering["first_n"].vars_kept)
        for name, counts in of_each.items()
    } == dict.fromkeys(of_each, (100, 100))


def test_filter_first_n_absent_leaves_stopped_early_false_in_every_consumer(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """With no filter of the first n every consumer says the pass was not
    ended, so the test above can fail."""
    of_each = _every_consumer(_many(reference_vcf_dir), tmp_path)

    assert {name: counts.stopped_early for name, counts in of_each.items()} == (
        dict.fromkeys(of_each, False)
    )
