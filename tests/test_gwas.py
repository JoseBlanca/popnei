"""The association study of a continuous trait, from Python.

`docs/specs/gwas.md` has what is computed, the `GWASResult` it comes in and
the numbers this file asserts. Every literal here is the spec's, or is read
from `tests/reference/gwas/`, which
`tests/reference/gwas/make_reference.py` writes: the trait, the two
covariates and the subpopulation of each of the 200 individuals in
`phenotypes.csv`, and what plink2 v2.0.0-a.7.7 answered for the panel with
every genotype called in `plink2.panel_called.glm.linear.tsv`, 1200 rows.
The panel itself is `tests/reference/kinship/panel_called.vcf.gz`.

The three checks of "How it is verified" of "What every model shares" that
are made at the Python `calc_gwas` are here: every column against plink2
over all 1200 variants, every column against pyNei at commit ef0ca6e, which
`pyproject.toml` names, and that the blocks the source gives change nothing.
The worked example and the six literals are cargo tests, where the spec puts
them. What the rest of this file covers is what only this layer has: which
individuals are tested and in which order, the design built out of a frame
of covariates, and what a wrong argument is refused with.

A `beta` and an `se` are compared with a share of the `se` of their variant
and never with a share of themselves, which is the rule of "How it is
verified" of "What every model shares": a study is mostly null, so for most
variants `beta` is what the rounding of a sum of cancelling products left,
and a bound relative to it asks for an accuracy no arithmetic has.
"""

import pathlib

import numpy
import pandas
import popnei
import pytest
from popnei import GWASModel, TraitType, calc_gwas, open_vars, open_vcf, write_vars
from pynei import vars_from_vcf
from pynei.gwas import calc_gwas as pynei_gwas

REFERENCE_GWAS_DIR = pathlib.Path(__file__).parent / "reference" / "gwas"
REFERENCE_KINSHIP_DIR = pathlib.Path(__file__).parent / "reference" / "kinship"

# The panel with every genotype called, which plink2 and pyNei were run on:
# 200 individuals and 1200 biallelic diploid variants on two chromosomes.
PANEL = REFERENCE_KINSHIP_DIR / "panel_called.vcf.gz"
PANEL_NUM_VARS = 1200
PANEL_NUM_INDIVIDUALS = 200

# How far a `beta` or an `se` of the panel may be from plink2's, as a share
# of the `se` plink2 printed for that variant. It is the 1e-5 of "How it is
# verified" of "The linear model" of the spec.
OF_PLINK2 = 1e-5

# How far a p-value of the panel may be from plink2's, as a share of it, and
# how far an `allele_freq` may be from plink2's `A1_FREQ`, absolute: a
# frequency lies between 0 and 1, so it is compared without a scale.
OF_PLINK2_P_VALUE = 1e-5
OF_PLINK2_FREQUENCY = 1e-6

# How many significant digits plink2 prints `BETA`, `SE` and `P` to, and so
# how much of the last of them is the most its printing rounds a value by:
# half a unit in that last place.
PLINK2_DIGITS = 6

# How far a `beta` or an `se` of the panel may be from pyNei's, as a share of
# the `se` pyNei gives for that variant, and how far a p-value may be, as the
# distance between the two in log10.
#
# "How it is verified" of "What every model shares" asks for 1e-9 relative
# and says that such a number is lowered until it fails and set two or three
# times above where it broke. Both were, on 24 September 2026 over all 1200
# variants on both backends. The `beta` breaks at 5.8e-15: the furthest from
# pyNei is `var0482`, 1.11e-15 away with an `se` of 0.190, which is
# 5.8296e-15 of it, the same variant and the same difference on Accelerate
# and on faer. So `OF_PYNEI` is 2.6 times the worst measured, where the
# spec's 1e-9 had 170000 times it. The `se` of every variant is nearer than
# that, 8.28e-16 of its own `se` at worst on Accelerate and 8.40e-16 on
# faer. The p-value breaks at 3.79e-12, at `var0216`, whose p-value of
# 0.99961 the two libraries give 8.74e-12 apart as a share of it, and
# `OF_PYNEI_P_VALUE` is 2.6 times that.
OF_PYNEI = 1.5e-14
OF_PYNEI_P_VALUE = 1e-11

# How far a column of the study of the panel read in blocks of 77 may be
# from the same column of the study of the panel read in one block, as a
# share of the value.
#
# The spec asks for 1e-12 and, as above, for the number to be lowered until
# it fails. This one fails at no number: the two studies give the same bits,
# a difference of 0 over all four columns and all 1200 variants on both
# backends on 24 September 2026. Nothing of the linear model is summed
# across the variants, so a variant's four numbers are the same whichever
# block it was read in, and the bound stays the spec's.
OF_THE_BLOCKS = 1e-12

# How many variants a batch of the vars file the panel is written to holds,
# which is what the pass over it reads at a time.
VARS_PER_BLOCK = 77

# The worked example of "The worked example" of the spec: six diploid
# individuals, one covariate beside the intercept and three variants, of
# which the third has every individual heterozygous and so no variance. The
# `./.` of `i3` at `v1` is the genotype that takes the mean dosage of its
# variant, 0.8.
WORKED_EXAMPLE_INDIVIDUALS = ("i0", "i1", "i2", "i3", "i4", "i5")
WORKED_EXAMPLE_GENOTYPES = (
    ("0/0", "0/1", "1/1", "0/0", "0/1", "1/1"),
    ("0/0", "0/1", "1/1", "./.", "0/1", "0/0"),
    ("0/1", "0/1", "0/1", "0/1", "0/1", "0/1"),
)
WORKED_EXAMPLE_TRAIT = (2.0, 3.0, 5.0, 4.0, 4.0, 7.0)
WORKED_EXAMPLE_COVARIATE = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

# Its null model, which pyNei gave at commit ef0ca6e: an intercept of
# 3.6666666666666683, an effect of the covariate of 1.0 and a residual
# variance of 3.333333333333334, which is its residual sum of squares of
# 13.333333333333336 over the 4 degrees of freedom of 6 individuals and 2
# coefficients.
WORKED_EXAMPLE_NULL = {
    "intercept": 3.6666666666666683,
    "cov": 1.0,
    "residual_variance": 3.333333333333334,
}

# Its three rows, which pyNei gave as well: `v2` keeps its frequency and has
# no test.
WORKED_EXAMPLE_ROWS = (
    ("v0", 0.5, 1.5, 0.600925212577332, 0.088004892382756),
    ("v1", 0.4, 0.3125, 1.305204592306424, 0.826200867452417),
    ("v2", 0.5, numpy.nan, numpy.nan, numpy.nan),
)

# How far a number of the worked example may be from pyNei's, as a share of
# that number: the cargo test of the same twelve numbers measures 1.24e-15
# on faer and 7.4e-16 on Accelerate and is set at 3e-15, and this one reads
# the same numbers through another layer. No number of the example is near
# 0, the smallest being the `beta` of `v1`, 0.3125 against an `se` of 1.3, so
# a share of the number and a share of the scale are the same bound here.
OF_THE_WORKED_EXAMPLE = 3e-15


def _phenotypes() -> pandas.DataFrame:
    """The trait, the covariates and the subpopulation of each individual of
    the panel, indexed by the name of the individual."""
    return pandas.read_csv(REFERENCE_GWAS_DIR / "phenotypes.csv", index_col=0)


def _of_plink2() -> pandas.DataFrame:
    """What plink2 `--glm hide-covar` answered for the panel, one row per
    variant in the order of the file."""
    return pandas.read_csv(
        REFERENCE_GWAS_DIR / "plink2.panel_called.glm.linear.tsv", sep="\t"
    )


def _the_study_of_the_panel(variants=None):
    """The study of the panel with `cov1` and `cov2`, which is what plink2
    and pyNei were given."""
    phenotypes = _phenotypes()
    return calc_gwas(
        open_vcf(PANEL) if variants is None else variants,
        phenotypes["cont"],
        TraitType.CONTINUOUS,
        covariates=phenotypes[["cov1", "cov2"]],
    )


def _half_of_the_last_digit(values: numpy.ndarray) -> numpy.ndarray:
    """The most that printing `values` to six significant digits rounds each
    of them by: half a unit in the place of the last digit printed.

    A value of 0.139354 is printed to the millionths, so its rounding is at
    most 5e-7, and one of 1.03892 to the hundred thousandths, 5e-6. That step
    at every power of ten is the whole of what is added to the tolerance
    below, and the reason it cannot be left out.
    """
    return 0.5 * 10.0 ** (
        numpy.floor(numpy.log10(numpy.abs(values))) - (PLINK2_DIGITS - 1)
    )


def _assert_within_the_printed_scale(
    ours: numpy.ndarray, theirs: numpy.ndarray, scale: numpy.ndarray, ids, what: str
) -> None:
    """Every value of `ours` within `OF_PLINK2` times `scale`, the `se` of
    its variant, plus half of the last digit plink2 printed for the value it
    is compared with.

    The second term is not slack: plink2 prints six significant digits, so
    what it wrote for a `beta` of 1.03892 is up to 5e-6 away from the number
    it computed, while `OF_PLINK2` times that variant's `se` of 0.190445 is
    1.9e-6. A tolerance is a budget shared between the rounding of the number
    it is compared against and the difference it is meant to catch, and here
    the first share is the larger of the two for the two variants of the
    panel whose `beta` passes 1.
    """
    budget = OF_PLINK2 * scale + _half_of_the_last_digit(theirs)
    difference = numpy.abs(ours - theirs)
    over = numpy.flatnonzero(difference > budget)
    assert over.size == 0, (
        f"{what} of {[ids[at] for at in over[:5]]} is further from plink2 than "
        f"the {OF_PLINK2} of the `se` plus half of the last digit printed: "
        f"the worst is {difference[over].max()} against a budget of "
        f"{budget[over[numpy.argmax(difference[over])]]}"
    )


def test_every_variant_of_the_panel_is_plink2s() -> None:
    """The four columns of all 1200 variants against plink2's.

    The ids are compared first, so that every row is matched to the variant
    plink2 wrote it for and not to the row at the same place.

    Measured over the 1200 on 24 September 2026, the same on Accelerate and
    on faer: every `allele_freq` is plink2's `A1_FREQ` exactly; the worst
    `beta` is 3.69e-6 away, at `var0482`, which is 53 per cent of the
    6.90e-6 of its budget, of which the `se` gives 1.90e-6 and the printed
    digits 5e-6; the worst `se` is 4.90e-7 away, at `var0680`, 30 per cent of
    its budget of 1.65e-6; and the worst p-value is 4.52e-6 of itself, 45 per
    cent of the 1e-5 allowed.

    As a share of the `se` alone the worst `beta` is 1.94e-5, which is why
    the printed digits are in the tolerance: `var0482` has a `beta` of
    1.03892, and six significant digits of a value above 1 are rounded by up
    to 5e-6 absolute, 2.6e-5 of its `se`. No arithmetic can meet a bound that
    the printing of the number it is compared against already passes.
    """
    result = _the_study_of_the_panel()
    of_plink2 = _of_plink2()

    assert len(result.stats.index) == PANEL_NUM_VARS
    assert len(of_plink2.index) == PANEL_NUM_VARS
    assert list(result.stats["id"]) == list(of_plink2["ID"])
    assert result.null_model.num_individuals == PANEL_NUM_INDIVIDUALS
    se = of_plink2["SE"].to_numpy()
    ids = list(of_plink2["ID"])
    numpy.testing.assert_allclose(
        result.stats["allele_freq"].to_numpy(),
        of_plink2["A1_FREQ"].to_numpy(),
        rtol=0,
        atol=OF_PLINK2_FREQUENCY,
    )
    _assert_within_the_printed_scale(
        result.stats["beta"].to_numpy(),
        of_plink2["BETA"].to_numpy(),
        se,
        ids,
        "the effect",
    )
    _assert_within_the_printed_scale(
        result.stats["se"].to_numpy(), se, se, ids, "the standard error"
    )
    numpy.testing.assert_allclose(
        result.stats["p_value"].to_numpy(),
        of_plink2["P"].to_numpy(),
        rtol=OF_PLINK2_P_VALUE,
        atol=0,
    )


def test_every_variant_of_the_panel_is_pyneis() -> None:
    """Both libraries on the same VCF, over all 1200 variants and the null
    model they were tested against.

    pyNei gives no counts of the pass, so what is compared is the four
    columns, the individuals that were tested, the coefficients of the null
    model and its residual variance. The variants that have no answer are
    asserted to be the same ones, which on this panel is none of them: every
    variant of it varies among the 200 individuals.
    """
    ours = _the_study_of_the_panel()
    phenotypes = _phenotypes()
    theirs = pynei_gwas(
        vars_from_vcf(PANEL),
        phenotypes["cont"],
        "continuous",
        covariates=phenotypes[["cov1", "cov2"]],
    )

    assert ours.individuals == tuple(theirs.samples)
    assert ours.null_model.model == theirs.null_model.model
    assert ours.trait == theirs.trait
    assert ours.test == theirs.test
    numpy.testing.assert_allclose(
        ours.null_model.covariate_effects.to_numpy(),
        theirs.null_model.covariate_effects.to_numpy(),
        rtol=OF_PYNEI,
        atol=0,
    )
    assert ours.null_model.residual_variance == pytest.approx(
        theirs.null_model.residual_variance, rel=OF_PYNEI
    )
    se = theirs.stats["se"].to_numpy()
    numpy.testing.assert_allclose(
        ours.stats["allele_freq"].to_numpy(),
        theirs.stats["allele_freq"].to_numpy(),
        rtol=0,
        atol=OF_PLINK2_FREQUENCY,
    )
    for column in ("beta", "se"):
        difference = numpy.abs(
            ours.stats[column].to_numpy() - theirs.stats[column].to_numpy()
        )
        assert numpy.nanmax(difference / se) <= OF_PYNEI, (
            f"the worst {column} is {numpy.nanmax(difference / se)} of the `se` "
            f"of its variant against the {OF_PYNEI} allowed"
        )
    ours_p = ours.stats["p_value"].to_numpy()
    theirs_p = theirs.stats["p_value"].to_numpy()
    assert (numpy.isnan(ours_p) == numpy.isnan(theirs_p)).all()
    numpy.testing.assert_allclose(
        numpy.log10(ours_p), numpy.log10(theirs_p), rtol=0, atol=OF_PYNEI_P_VALUE
    )


def test_the_blocks_the_source_gives_change_nothing(tmp_path: pathlib.Path) -> None:
    """The panel read in blocks of 77 gives the study of the panel read in
    one block.

    The VCF of 1200 variants is one block, since popnei reads 10000 of 200
    individuals at a time, and the same variants written to a vars file with
    77 variants a batch are 16 batches, 15 of 77 and one of 45, which is what
    the pass over that file reads at a time: the reader of a vars file builds
    a batch whole whatever size its blocks are asked for. Read with pyarrow
    on 24 September 2026, and it cannot be read through popnei, since every
    pass puts a `reblock` over its reader and gives the blocks that one
    makes.

    What the pass does with those 16 is join them into the block the study
    reads, so what this says is that a source which gives its variants a few
    at a time is studied as one that gives them all at once. That the study's
    own loop over its blocks adds the rows of each one after the rows of the
    ones before, and that a variant has the same answer in the second block
    as in the first, is the cargo test of 10100 variants, which is the size
    at which that loop runs twice.

    Measured on 24 September 2026 on both backends: the two agree to the bit
    over all four columns and all 1200 variants, a difference of 0.
    """
    of_the_vcf = _the_study_of_the_panel()
    path = tmp_path / "panel.vars"
    written = write_vars(open_vcf(PANEL), path, num_vars_per_block=VARS_PER_BLOCK)
    of_the_blocks = _the_study_of_the_panel(open_vars(path))

    assert written.pass_stats.num_vars == PANEL_NUM_VARS
    assert of_the_blocks.pass_stats.num_vars == PANEL_NUM_VARS
    assert list(of_the_blocks.stats["id"]) == list(of_the_vcf.stats["id"])
    for column in ("allele_freq", "beta", "se", "p_value"):
        numpy.testing.assert_allclose(
            of_the_blocks.stats[column].to_numpy(),
            of_the_vcf.stats[column].to_numpy(),
            rtol=OF_THE_BLOCKS,
            atol=0,
        )


def _worked_example_vcf(path: pathlib.Path, individuals=WORKED_EXAMPLE_INDIVIDUALS):
    """The worked example of the spec as a VCF at `path`: three variants of
    six diploid individuals, each variant with its own id and position."""
    lines = [
        "##fileformat=VCFv4.4",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + list(individuals)
        ),
    ]
    for at, genotypes in enumerate(WORKED_EXAMPLE_GENOTYPES):
        lines.append(
            "\t".join(
                [
                    "chr1",
                    str((at + 1) * 1000),
                    f"v{at}",
                    "A",
                    "T",
                    ".",
                    "PASS",
                    ".",
                    "GT",
                ]
                + list(genotypes)
            )
        )
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def worked_example(tmp_path: pathlib.Path) -> pathlib.Path:
    """The VCF of the worked example of the spec."""
    return _worked_example_vcf(tmp_path / "worked_example.vcf")


def _the_trait(individuals=WORKED_EXAMPLE_INDIVIDUALS) -> pandas.Series:
    """The trait of the worked example, indexed by individual."""
    return pandas.Series(WORKED_EXAMPLE_TRAIT, index=list(individuals))


def _the_covariate(individuals=WORKED_EXAMPLE_INDIVIDUALS) -> pandas.DataFrame:
    """Its one covariate beside the intercept, in a frame of one column."""
    return pandas.DataFrame({"cov": WORKED_EXAMPLE_COVARIATE}, index=list(individuals))


def _the_worked_example(path: pathlib.Path, **asked_for):
    """The study of the worked example, with `asked_for` over its trait and
    its covariate."""
    arguments = {
        "phenotype": _the_trait(),
        "trait": "continuous",
        "covariates": _the_covariate(),
    }
    arguments.update(asked_for)
    return calc_gwas(open_vcf(path, only_passed=False), **arguments)


def test_the_worked_example_gives_pyneis_null_model_and_its_three_rows(
    worked_example: pathlib.Path,
) -> None:
    """The result objects of the study, against the numbers pyNei gave.

    The calculation is the core's and is asserted there. What this adds is
    that they reach Python whole: the effects under the names of the columns
    of the design, the three columns that name a variant beside the four that
    answer for it, and the variant that has no variance with its frequency
    and three NaNs.
    """
    result = _the_worked_example(worked_example)

    assert result.null_model.model == GWASModel.LM
    assert result.trait == TraitType.CONTINUOUS
    assert result.test == popnei.TestType.WALD
    assert result.used_grammar_gamma_approx is False
    assert result.null_model.num_individuals == 6
    assert result.individuals == WORKED_EXAMPLE_INDIVIDUALS
    assert result.null_model.genetic_variance is None
    assert result.null_model.heritability is None
    assert list(result.null_model.covariate_effects.index) == ["intercept", "cov"]
    for name in ("intercept", "cov"):
        assert result.null_model.covariate_effects[name] == pytest.approx(
            WORKED_EXAMPLE_NULL[name], rel=OF_THE_WORKED_EXAMPLE
        )
    assert result.null_model.residual_variance == pytest.approx(
        WORKED_EXAMPLE_NULL["residual_variance"], rel=OF_THE_WORKED_EXAMPLE
    )
    assert list(result.stats.columns) == [
        "chrom",
        "pos",
        "id",
        "allele_freq",
        "beta",
        "se",
        "p_value",
    ]
    assert list(result.stats["chrom"]) == ["chr1"] * 3
    assert list(result.stats["pos"]) == [1000, 2000, 3000]
    assert list(result.stats["id"]) == ["v0", "v1", "v2"]
    for at, (var, allele_freq, beta, se, p_value) in enumerate(WORKED_EXAMPLE_ROWS):
        row = result.stats.iloc[at]
        assert row["id"] == var
        assert row["allele_freq"] == pytest.approx(
            allele_freq, rel=OF_THE_WORKED_EXAMPLE
        )
        if numpy.isnan(beta):
            assert numpy.isnan(row["beta"])
            assert numpy.isnan(row["se"])
            assert numpy.isnan(row["p_value"])
            continue
        assert row["beta"] == pytest.approx(beta, rel=OF_THE_WORKED_EXAMPLE)
        assert row["se"] == pytest.approx(se, rel=OF_THE_WORKED_EXAMPLE)
        assert row["p_value"] == pytest.approx(p_value, rel=OF_THE_WORKED_EXAMPLE)


def test_an_individual_with_no_phenotype_is_not_tested(
    worked_example: pathlib.Path,
) -> None:
    """The individuals tested are those with a phenotype, and the frequency
    of a variant is of them and not of the panel.

    The dosages of `v0` over the five that are left are 0, 1, 2, 0 and 1,
    whose mean is 0.8, so its `allele_freq` is 0.4 where the six give 0.5.
    """
    trait = _the_trait()
    trait["i5"] = numpy.nan

    result = _the_worked_example(worked_example, phenotype=trait)

    assert result.individuals == ("i0", "i1", "i2", "i3", "i4")
    assert result.null_model.num_individuals == 5
    assert result.stats["allele_freq"][0] == pytest.approx(
        0.4, rel=OF_THE_WORKED_EXAMPLE
    )


def test_the_individuals_are_tested_in_the_order_the_source_has_them(
    worked_example: pathlib.Path,
) -> None:
    """A phenotype written in another order gives the same study.

    The phenotype, the rows of the design and the dosages of a variant are
    read together row by row, so a study that took them in the order the
    phenotype was written in would measure one individual's trait against
    another's genotypes. Every number here would move if it did: the trait of
    the example rises with the dosages of `v0`.
    """
    turned = _the_trait().iloc[::-1]
    covariates = _the_covariate().iloc[::-1]

    result = _the_worked_example(
        worked_example, phenotype=turned, covariates=covariates
    )
    as_they_come = _the_worked_example(worked_example)

    assert result.individuals == WORKED_EXAMPLE_INDIVIDUALS
    for column in ("allele_freq", "beta", "se", "p_value"):
        numpy.testing.assert_allclose(
            result.stats[column].to_numpy(),
            as_they_come.stats[column].to_numpy(),
            rtol=0,
            atol=0,
        )


def test_the_study_is_over_the_variants_the_steps_keep(
    worked_example: pathlib.Path,
) -> None:
    """A filter on the `Variants` takes variants out of the study, and the
    counts of the pass say what it was given and kept.

    `v2` is the variant where every individual is heterozygous, which the
    filter by observed heterozygosity is the one that catches: it has no
    variance and no answer, and a study that keeps it reports its frequency
    alone.
    """
    variants = open_vcf(worked_example, only_passed=False)
    variants.filter_by_obs_het(max_allowed_obs_het=0.9)

    result = _the_worked_example(worked_example, phenotype=_the_trait())
    filtered = calc_gwas(
        variants,
        _the_trait(),
        "continuous",
        covariates=_the_covariate(),
    )

    assert list(result.stats["id"]) == ["v0", "v1", "v2"]
    assert list(filtered.stats["id"]) == ["v0", "v1"]
    assert filtered.pass_stats.num_vars == 2
    assert filtered.pass_stats.filtering["obs_het"].vars_processed == 3
    assert filtered.pass_stats.filtering["obs_het"].vars_kept == 2


def test_a_phenotype_of_a_frame_of_one_column_is_that_column(
    worked_example: pathlib.Path,
) -> None:
    """`frame[["cont"]]` is taken as the trait it holds, and a frame of more
    columns is refused as naming no one trait."""
    one_column = pandas.DataFrame(
        {"cont": WORKED_EXAMPLE_TRAIT}, index=list(WORKED_EXAMPLE_INDIVIDUALS)
    )

    of_the_frame = _the_worked_example(worked_example, phenotype=one_column)
    of_the_series = _the_worked_example(worked_example)

    numpy.testing.assert_allclose(
        of_the_frame.stats["beta"].to_numpy(),
        of_the_series.stats["beta"].to_numpy(),
        rtol=0,
        atol=0,
    )
    with pytest.raises(ValueError, match="frame of 2 columns"):
        _the_worked_example(
            worked_example, phenotype=one_column.assign(other=WORKED_EXAMPLE_TRAIT)
        )


def test_a_binomial_trait_is_refused_with_the_model_that_is_being_written(
    worked_example: pathlib.Path,
) -> None:
    """The logistic regression is not written, and the study that needs it
    says so."""
    binomial = pandas.Series(
        [0.0, 1.0, 0.0, 1.0, 0.0, 1.0], index=list(WORKED_EXAMPLE_INDIVIDUALS)
    )

    with pytest.raises(ValueError, match="logistic regression, which is being written"):
        _the_worked_example(worked_example, phenotype=binomial, trait="binomial")


def test_a_trait_of_another_name_is_refused_with_the_two_names(
    worked_example: pathlib.Path,
) -> None:
    """A name that is of neither trait is refused with both of them, which
    the binding crate holds because the core has them."""
    with pytest.raises(ValueError, match="`continuous`.*`binomial`.*`quantitative`"):
        _the_worked_example(worked_example, trait="quantitative")


def test_a_trait_that_is_no_name_is_refused_by_its_type(
    worked_example: pathlib.Path,
) -> None:
    """What a user writes for `trait` instead of a name is usually the
    covariates, which the argument before it takes."""
    with pytest.raises(TypeError, match="`trait` is"):
        _the_worked_example(worked_example, trait=_the_covariate())


@pytest.mark.parametrize(
    ("argument", "value"),
    [
        ("kinship", "a matrix"),
        ("test", "wald"),
        ("use_grammar_gamma_approx", True),
    ],
)
def test_what_the_linear_mixed_model_brings_is_refused_by_name(
    worked_example: pathlib.Path, argument: str, value: object
) -> None:
    """A kinship, a test and the approximation are refused by name.

    A call that gave a kinship and got a study without one would be a linear
    model reported as a mixed one, with nothing to show it, so they are
    refused and not ignored while the linear mixed model is being written.
    """
    with pytest.raises(
        ValueError, match=rf"`{argument}` belongs to the linear mixed model"
    ):
        _the_worked_example(worked_example, **{argument: value})


def test_an_individual_of_the_phenotype_that_the_dataset_has_not_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """A name of the phenotype that is of nobody the pass gives is refused by
    that name: a phenotype and a VCF that were matched by hand, and one name
    spelled wrong, would otherwise give a study of fewer individuals."""
    trait = _the_trait()
    trait["i9"] = 5.0

    with pytest.raises(ValueError, match="'i9' has a phenotype"):
        _the_worked_example(worked_example, phenotype=trait)


def test_an_individual_named_twice_in_the_phenotype_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """One individual with two phenotypes is refused, naming the two places
    it is at: it would weigh twice in the null model and in every variant."""
    trait = _the_trait()
    twice = pandas.concat([trait, trait.iloc[:1]])

    with pytest.raises(ValueError, match="'i0' is named twice"):
        _the_worked_example(worked_example, phenotype=twice)


def test_a_phenotype_that_is_not_a_number_is_refused_by_the_individual(
    worked_example: pathlib.Path,
) -> None:
    """A trait whose values are names is not a trait of a study, and the one
    that is not a number is named with its individual."""
    trait = _the_trait().astype(object)
    trait["i2"] = "tall"

    with pytest.raises(ValueError, match="the phenotype of 'i2' is 'tall'"):
        _the_worked_example(worked_example, phenotype=trait)


def test_a_covariate_that_does_not_cover_a_tested_individual_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """A covariate holds a value for every individual that is tested."""
    covariates = _the_covariate().drop(index="i5")

    with pytest.raises(ValueError, match="the covariate 'cov' has no value for 'i5'"):
        _the_worked_example(worked_example, covariates=covariates)


def test_a_covariate_with_a_missing_value_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """A covariate that is there and empty is refused where one that is not
    there at all is, and the message says what to do with it."""
    covariates = _the_covariate()
    covariates.loc["i3", "cov"] = numpy.nan

    with pytest.raises(ValueError, match="at 'i3' is missing"):
        _the_worked_example(worked_example, covariates=covariates)


def test_a_covariate_whose_values_are_names_is_refused_with_what_to_write(
    worked_example: pathlib.Path,
) -> None:
    """A categorical covariate is refused with the one column per value that
    `pandas.get_dummies` writes, which is what a user does with it."""
    covariates = _the_covariate().astype(object)
    covariates["cov"] = ["north", "south"] * 3

    with pytest.raises(ValueError, match="pandas.get_dummies"):
        _the_worked_example(worked_example, covariates=covariates)


def test_a_covariate_named_intercept_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """`intercept` is the name the effect of the column of ones comes back
    under, so a covariate of that name would be one row of
    `covariate_effects` with it and the user would read one of the two
    without knowing which."""
    covariates = _the_covariate().rename(columns={"cov": "intercept"})

    with pytest.raises(ValueError, match="a covariate is named 'intercept'"):
        _the_worked_example(worked_example, covariates=covariates)


def test_a_covariate_that_is_a_copy_of_another_is_refused_as_collinear(
    worked_example: pathlib.Path,
) -> None:
    """Covariates that are not independent have no one set of effects, and
    the fit would answer with whichever the arithmetic reached."""
    covariates = _the_covariate()
    covariates["twice"] = covariates["cov"]

    with pytest.raises(
        ValueError, match="3 columns, the intercept among them, and only 2"
    ):
        _the_worked_example(worked_example, covariates=covariates)


def test_a_covariate_that_is_not_finite_is_refused_by_its_place(
    worked_example: pathlib.Path,
) -> None:
    """An infinity is a number to pandas and not to a fit, and the core names
    the column of the design and the individual it is in.

    It is what a covariate that came out of the user's own arithmetic as an
    infinity gives, a division by 0 among the causes, and it names no file:
    what is wrong is wrong whatever variants are read.
    """
    covariates = _the_covariate()
    covariates.loc["i1", "cov"] = numpy.inf

    with pytest.raises(ValueError, match="tested individual 1 is inf") as refused:
        _the_worked_example(worked_example, covariates=covariates)
    assert str(refused.value).startswith("the value of the column 1")


def test_a_study_of_too_few_individuals_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """A study needs the columns of its design plus two individuals, or the
    variant leaves nothing to measure its uncertainty from."""
    trait = _the_trait()
    trait[["i3", "i4", "i5"]] = numpy.nan

    with pytest.raises(
        ValueError, match="3 individuals are tested and the design has 2 columns"
    ):
        _the_worked_example(worked_example, phenotype=trait)


def test_what_is_not_a_variants_is_refused_by_its_type() -> None:
    """The path of the VCF is what a user gives instead, and what that gave
    before this check was the `AttributeError` of an object with no source
    inside it."""
    with pytest.raises(TypeError, match="`calc_gwas` reads the variants"):
        calc_gwas(str(PANEL), _the_trait(), "continuous")
