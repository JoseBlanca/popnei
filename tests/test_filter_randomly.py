"""The filter that keeps variants at random from Python: which variants it
keeps, that every pass of one `Variants` keeps the same ones, and what it
refuses at the call.

`docs/specs/filters.md` has the filter under "The filter that keeps variants
at random". pyNei has no such filter, so the variants are those of
`tests/reference/filters/random_draws.py`, a Python version of the rule
whose draws were checked against those of Java's
`java.util.SplittableRandom`, applied to the 500 variants of `many.vcf` of
`docs/specs/io_vcf.md` read with every variant given, those that failed their
FILTER among them, and run on 5 October 2026.
"""

import math
from pathlib import Path

import pandas
import pytest
from popnei import (
    FilteringStats,
    PassStats,
    Step,
    calc_gwas,
    calc_kinship,
    calc_pairwise_kosman_dists,
    do_pca_from_variants,
    open_vars,
    open_vcf,
    write_vars,
)
from popnei.variant import Variants

MANY_NUM_VARS = 500

# The three rows of the table of the spec: a keep rate and a seed, how many
# of the 500 variants they keep and the first five kept, by position on chr1.
KEPT_AT_0_1_AND_42 = 45
FIRST_FIVE_AT_0_1_AND_42 = [1148, 1666, 1777, 1888, 2332]
TABLE_OF_THE_SPEC = [
    (0.1, 42, KEPT_AT_0_1_AND_42, FIRST_FIVE_AT_0_1_AND_42),
    (0.5, 42, 243, [1037, 1074, 1111, 1148, 1222]),
    (0.1, 7, 49, [1037, 1962, 2147, 2332, 2591]),
]


def _many(reference_vcf_dir: Path) -> Variants:
    """The 500 variants of `many.vcf`, the ones that failed their FILTER
    among them, which are those the numbers of the spec are of."""
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


def _counts_of_the_45() -> PassStats:
    """The counts of a pass over `many.vcf` with the filter at 0.1 and a seed
    of 42."""
    return PassStats(
        num_vars=KEPT_AT_0_1_AND_42,
        filtering={
            "random": FilteringStats(
                vars_processed=MANY_NUM_VARS, vars_kept=KEPT_AT_0_1_AND_42
            )
        },
        stopped_early=False,
    )


@pytest.mark.parametrize(
    ("keep_rate", "seed", "num_kept", "first_five"), TABLE_OF_THE_SPEC
)
@pytest.mark.parametrize("num_vars_per_block", [7, None])
def test_filter_randomly_keeps_the_variants_of_the_table_of_the_spec(
    reference_vcf_dir: Path,
    keep_rate: float,
    seed: int,
    num_vars_per_block: int | None,
    num_kept: int,
    first_five: list[int],
) -> None:
    """Each row in blocks of 7 and of the default size: the draws are one
    for each variant, so the size of the blocks changes nothing, and the
    third row is there so that a seed that is not read fails."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(keep_rate, seed)

    kept, pass_stats = _kept(variants, num_vars_per_block)

    assert len(kept) == num_kept
    assert kept[:5] == [("chr1", pos) for pos in first_five]
    assert pass_stats.filtering["random"] == FilteringStats(
        vars_processed=MANY_NUM_VARS, vars_kept=num_kept
    )


def test_filter_randomly_gives_the_same_45_in_two_iter_blocks(
    reference_vcf_dir: Path,
) -> None:
    """Every pass starts the generator at the seed again, so a second pass
    keeps the variants the first one kept."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1, 42)

    first, of_the_first = _kept(variants)
    second, of_the_second = _kept(variants, 7)

    assert first[:5] == [("chr1", pos) for pos in FIRST_FIVE_AT_0_1_AND_42]
    assert second == first
    assert of_the_first == of_the_second == _counts_of_the_45()


def test_filter_randomly_gives_the_45_to_the_distances_and_to_the_pca(
    reference_vcf_dir: Path,
) -> None:
    """Two calculations on one `Variants` see one sample."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1)

    kosman = calc_pairwise_kosman_dists(variants)
    pca = do_pca_from_variants(variants, transform_to_biallelic=True)

    assert kosman.pass_stats == _counts_of_the_45()
    assert pca.pass_stats == _counts_of_the_45()


def _the_45_in_a_vars_file(reference_vcf_dir: Path, tmp_path: Path) -> Variants:
    """The 45 variants written to a vars file, opened with no step."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1)
    path = tmp_path / "the_45.vars"
    written = write_vars(variants, path)
    assert written.pass_stats == _counts_of_the_45()
    return open_vars(path)


def test_filter_randomly_gives_a_pca_whose_second_pass_reads_the_same_45(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """With 10 components the PCA reads the source a second time for the
    weights of the variants. That pass does not fail with the error of a
    second pass that gave other variants, and the result is the one of the
    45 written to a file with no step: a second pass that read other
    variants would give other weights. The two analyses are given the same
    dosages in the same order, so the numbers are compared exactly.
    `many.vcf` holds variants of three alleles, so the dosages count every
    allele that is not the major one."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1)
    over_the_file = _the_45_in_a_vars_file(reference_vcf_dir, tmp_path)

    of_the_filter = do_pca_from_variants(
        variants, transform_to_biallelic=True, num_prin_comps=10
    )
    of_the_file = do_pca_from_variants(
        over_the_file, transform_to_biallelic=True, num_prin_comps=10
    )

    assert of_the_filter.pass_stats == _counts_of_the_45()
    assert of_the_filter.princomps.shape == (10, KEPT_AT_0_1_AND_42)
    pandas.testing.assert_frame_equal(
        of_the_filter.projections, of_the_file.projections, check_exact=True
    )
    pandas.testing.assert_frame_equal(
        of_the_filter.princomps, of_the_file.princomps, check_exact=True
    )
    pandas.testing.assert_series_equal(
        of_the_filter.explained_variance_percent,
        of_the_file.explained_variance_percent,
        check_exact=True,
    )


def test_filter_randomly_gives_a_gwas_whose_grammar_gamma_pass_reads_the_same_45(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The GRAMMAR-Gamma approximation reads the source a second time with no
    check of its own, so the study over the filter equals the one over the
    45 written to a file only when both of its passes saw the 45. The two are
    given the same dosages, phenotype and kinship, so they are compared
    exactly."""
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1)
    over_the_file = _the_45_in_a_vars_file(reference_vcf_dir, tmp_path)
    kinship = calc_kinship(_many(reference_vcf_dir), transform_to_biallelic=True)
    individuals = variants.individuals
    phenotype = pandas.Series(
        [float(index % 7) for index in range(len(individuals))],
        index=list(individuals),
    )

    of_the_filter = calc_gwas(
        variants,
        phenotype,
        "continuous",
        kinship=kinship,
        use_grammar_gamma_approx=True,
        transform_to_biallelic=True,
    )
    of_the_file = calc_gwas(
        over_the_file,
        phenotype,
        "continuous",
        kinship=kinship,
        use_grammar_gamma_approx=True,
        transform_to_biallelic=True,
    )

    assert of_the_filter.used_grammar_gamma_approx is True
    assert of_the_filter.pass_stats == _counts_of_the_45()
    assert len(of_the_filter.stats) == KEPT_AT_0_1_AND_42
    pandas.testing.assert_frame_equal(
        of_the_filter.stats, of_the_file.stats, check_exact=True
    )
    assert of_the_filter.stats["p_value"].notna().all()


@pytest.mark.parametrize(
    ("keep_rate", "written"), [(-0.1, "-0.1"), (1.5, "1.5"), (math.nan, "NaN")]
)
def test_filter_randomly_refuses_a_keep_rate_out_of_0_to_1(
    reference_vcf_dir: Path, keep_rate: float, written: str
) -> None:
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError, match=rf"`keep_rate` is {written}, .* from 0 to 1"):
        variants.filter_randomly(keep_rate)
    assert variants.steps == ()


@pytest.mark.parametrize("keep_rate", [True, "0.1", None])
def test_filter_randomly_refuses_a_keep_rate_that_is_no_number_with_a_type_error(
    reference_vcf_dir: Path, keep_rate: object
) -> None:
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError, match="`keep_rate`"):
        variants.filter_randomly(keep_rate)
    assert variants.steps == ()


@pytest.mark.parametrize("seed", [-1, 2**64])
def test_filter_randomly_refuses_a_seed_that_64_bits_do_not_hold(
    reference_vcf_dir: Path, seed: int
) -> None:
    """A negative seed is the `ValueError` that names `seed`, and not the
    `OverflowError` of the conversion."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(ValueError, match=f"`seed` is {seed}, .* 18446744073709551615"):
        variants.filter_randomly(0.1, seed)
    assert variants.steps == ()


@pytest.mark.parametrize("seed", [1.5, True, "42", None])
def test_filter_randomly_refuses_a_seed_that_is_no_whole_number_with_a_type_error(
    reference_vcf_dir: Path, seed: object
) -> None:
    """`True` would be the seed 1 with nothing said."""
    variants = _many(reference_vcf_dir)

    with pytest.raises(TypeError, match="`seed`"):
        variants.filter_randomly(0.1, seed)
    assert variants.steps == ()


def test_filter_randomly_takes_the_largest_seed_of_64_bits(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1, 2**64 - 1)

    assert variants.steps == (
        Step(kind="random", args={"keep_rate": 0.1, "seed": 2**64 - 1}),
    )


def test_filter_randomly_refuses_a_second_filter_of_its_kind(
    reference_vcf_dir: Path,
) -> None:
    variants = _many(reference_vcf_dir)
    variants.filter_randomly(0.1)

    with pytest.raises(
        ValueError,
        match="by random already, with a keep rate of 0.1 and a seed of 42, and a "
        "second filter of that kind, with a keep rate of 0.2 and a seed of 7,",
    ):
        variants.filter_randomly(0.2, 7)
    assert variants.steps == (Step(kind="random", args={"keep_rate": 0.1, "seed": 42}),)


def test_filter_randomly_is_refused_after_the_filter_of_the_first_n(
    reference_vcf_dir: Path,
) -> None:
    """It would leave fewer than the n that filter keeps."""
    variants = _many(reference_vcf_dir)
    variants.filter_first_n(10)

    with pytest.raises(ValueError) as refusal:
        variants.filter_randomly(0.1)

    assert "a filter by random after it" in str(refusal.value)
    assert "first_n already" in str(refusal.value)
    assert [step.kind for step in variants.steps] == ["first_n"]


def test_filter_randomly_is_a_step_of_its_kind_and_both_arguments(
    reference_vcf_dir: Path,
) -> None:
    """Without a seed the step has the default, 42; a keep rate written as a
    whole number is read back as a float."""
    with_the_default = _many(reference_vcf_dir)
    with_the_default.filter_randomly(0.1)
    with_a_seed = _many(reference_vcf_dir)
    with_a_seed.filter_randomly(1, seed=7)

    assert with_the_default.steps == (
        Step(kind="random", args={"keep_rate": 0.1, "seed": 42}),
    )
    assert with_a_seed.steps == (
        Step(kind="random", args={"keep_rate": 1.0, "seed": 7}),
    )
    assert isinstance(with_a_seed.steps[0].args["keep_rate"], float)
    assert "random(keep_rate=0.1, seed=42)" in repr(with_the_default)
