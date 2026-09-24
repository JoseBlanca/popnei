"""The association study of a continuous trait, from Python.

`docs/specs/gwas.md` has what is computed, the `GWASResult` it comes in and
the numbers this file asserts. Every literal here is the spec's, or is read
from `tests/reference/gwas/`, which
`tests/reference/gwas/make_reference.py` writes: the trait, the two
covariates and the subpopulation of each of the 200 individuals in
`phenotypes.csv`, what plink2 v2.0.0-a.7.7 answered for the panel with every
genotype called in `plink2.panel_called.glm.linear.tsv`, and what GMMAT
1.5.0 and rrBLUP 4.6.3 answered for the mixed model in
`gmmat.null_models.tsv`, `gmmat.panel_called.lmm.score.tsv`,
`gmmat.panel.lmm.score.tsv` and `rrblup.panel_called.lmm.tsv`, 1200 rows
each. The panel itself is `tests/reference/kinship/panel_called.vcf.gz`, and
the kinship the mixed model is given is the one plink2 wrote for it, which
came from neither popnei nor pyNei.

The checks of "How it is verified" of "What every model shares" that are
made at the Python `calc_gwas` are here: every column against the reference
program of each model over all 1200 variants, every column against pyNei at
commit ef0ca6e, which `pyproject.toml` names, that the blocks the source
gives change nothing, and that the mixed model finds the variants that were
planted. The worked example and the six literals of each model are cargo
tests, where the spec puts them. What the rest of this file covers is what
only this layer has: which individuals are tested and in which order, the
design built out of a frame of covariates, the kinship cut to the tested
individuals, and what a wrong argument is refused with.

A `beta` and an `se` are compared with a share of the `se` of their variant
and never with a share of themselves, which is the rule of "How it is
verified" of "What every model shares": a study is mostly null, so for most
variants `beta` is what the rounding of a sum of cancelling products left,
and a bound relative to it asks for an accuracy no arithmetic has.
"""

import gzip
import json
import pathlib

import numpy
import pandas
import popnei
import pytest
from popnei import (
    GWASModel,
    Kinship,
    TraitType,
    _core,
    calc_gwas,
    open_vars,
    open_vcf,
    write_vars,
)
from pynei import vars_from_vcf
from pynei.gwas import Kinship as PyneiKinship
from pynei.gwas import calc_gwas as pynei_gwas

REFERENCE_GWAS_DIR = pathlib.Path(__file__).parent / "reference" / "gwas"
REFERENCE_KINSHIP_DIR = pathlib.Path(__file__).parent / "reference" / "kinship"

# The panel with every genotype called, which plink2 and pyNei were run on:
# 200 individuals and 1200 biallelic diploid variants on two chromosomes.
PANEL = REFERENCE_KINSHIP_DIR / "panel_called.vcf.gz"
PANEL_NUM_VARS = 1200
PANEL_NUM_INDIVIDUALS = 200

# The same 200 individuals and 1200 variants with 3 in 100 of the genotypes
# missing whole, which is the panel of `docs/specs/dists.md` and the one
# where a genotype takes the mean dosage of its variant. plink2 was not run
# on it, so it is compared with pyNei alone.
PANEL_WITH_MISSING = (
    pathlib.Path(__file__).parent / "reference" / "dists" / "panel.vcf.gz"
)
PANELS = (PANEL, PANEL_WITH_MISSING)

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

# How far a `beta` or an `se` of a panel may be from pyNei's, as a share of
# the `se` pyNei gives for that variant, and how far a p-value may be, as the
# distance between the two in log10.
#
# "How it is verified" of "What every model shares" asks for 1e-9 relative
# and says that such a number is lowered until it fails and set two or three
# times above where it broke. Both were, on 24 September 2026 over the 1200
# variants of each panel on both backends. The `beta` breaks at 6.8e-15: the
# furthest from pyNei is `var0671` of the panel with genotypes missing, on
# faer, 6.796e-15 of its `se`, and the worst on Accelerate is `var0482` of
# `panel_called`, 1.11e-15 away with an `se` of 0.190, 5.830e-15 of it. So
# `OF_PYNEI` is 2.2 times the worst measured, where the spec's 1e-9 had
# 150000 times it. Every `se` is nearer than that, 1.24e-15 of its own `se`
# at worst. The p-value breaks at 3.94e-12, at `var0717` of the panel with
# genotypes missing, and `OF_PYNEI_P_VALUE` is 2.5 times that.
OF_PYNEI = 1.5e-14
OF_PYNEI_P_VALUE = 1e-11

# How many variants a batch of the vars file the panel is written to holds,
# which is what the reader of that file gives the pass at a time.
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


@pytest.mark.parametrize("panel", PANELS)
def test_every_variant_of_a_panel_is_pyneis(panel: pathlib.Path) -> None:
    """Both libraries on the same VCF, over all 1200 variants of each panel
    and the null model they were tested against.

    pyNei gives no counts of the pass, so what is compared is the four
    columns, the individuals that were tested, the coefficients of the null
    model and its residual variance. The variants that have no answer are
    asserted to be the same ones, which on both panels is none of them.

    Both panels are here because the rule that gives a genotype with any
    allele missing the mean dosage of its variant runs at scale on one of
    them and on no other test of this file: `panel_called` has every
    genotype called, and `panel` has 3 in 100 of them missing whole.
    Measured on 24 September 2026, the worst `beta` of `panel_called` is
    `var0482` at 5.830e-15 of its `se` on both backends, and of `panel` it
    is `var0482` at 5.804e-15 on Accelerate and `var0671` at 6.796e-15 on
    faer. Neither panel has a variant that popnei or pyNei leaves without an
    answer, so the two sets of NaN are both empty and the check that they
    are the same ones is the weakest of the three; the five individuals of
    `test_a_frequency_of_half_called_genotypes_can_pass_a_half` are where a
    variant with no answer is asserted.
    """
    ours = _the_study_of_the_panel(open_vcf(panel))
    phenotypes = _phenotypes()
    theirs = pynei_gwas(
        vars_from_vcf(panel),
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


def test_the_vars_reader_gives_the_study_the_vcf_reader_gives(
    tmp_path: pathlib.Path,
) -> None:
    """The panel written to a vars file of 77 variants a batch is the study
    the VCF gives, to the bit.

    This is not a test of the size of the blocks the study reads, although
    it was written as one. Every pass puts a `reblock` over its reader, so
    the 16 batches of the file, 15 of 77 and one of 45, are joined into the
    one block of 1200 that the study reads, and the VCF gives that block
    whole; the two runs are the same arithmetic over the same values, and
    anything else would be a defect of one of the two readers rather than a
    digit moved. What runs the study's loop over more than one block is the
    cargo test of 10100 variants, which is where the size of a block is
    asserted.

    So what is compared is the bits, with no tolerance: a difference of 0
    over the four columns and all 1200 variants, measured on both backends
    on 24 September 2026. A tolerance here would pass a study that had read
    the file wrongly.
    """
    of_the_vcf = _the_study_of_the_panel()
    path = tmp_path / "panel.vars"
    written = write_vars(open_vcf(PANEL), path, num_vars_per_block=VARS_PER_BLOCK)
    of_the_vars_file = _the_study_of_the_panel(open_vars(path))

    assert written.pass_stats.num_vars == PANEL_NUM_VARS
    assert of_the_vars_file.pass_stats.num_vars == PANEL_NUM_VARS
    assert list(of_the_vars_file.stats["id"]) == list(of_the_vcf.stats["id"])
    assert list(of_the_vars_file.stats["chrom"]) == list(of_the_vcf.stats["chrom"])
    assert list(of_the_vars_file.stats["pos"]) == list(of_the_vcf.stats["pos"])
    for column in ("allele_freq", "beta", "se", "p_value"):
        ours = of_the_vars_file.stats[column].to_numpy()
        theirs = of_the_vcf.stats[column].to_numpy()
        assert (ours == theirs).all(), (
            f"{column} is not the same bits read from a vars file of "
            f"{VARS_PER_BLOCK} variants a batch as read from the VCF: the "
            f"worst is {numpy.abs(ours - theirs).max()}"
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


def test_the_approximation_of_a_study_with_no_kinship_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """The GRAMMAR-Gamma approximation stands in for the denominator of a
    mixed model's test, and a study with no kinship has no such denominator
    to approximate.

    It is the core that refuses it and not this layer, so the message is the
    one TypeScript gives, and it is not the message of a study that brought a
    kinship: there the denominator exists and popnei cannot approximate it
    yet, which the next test asserts.
    """
    with pytest.raises(ValueError, match="has no such denominator"):
        _the_worked_example(worked_example, use_grammar_gamma_approx=True)


def test_the_approximation_of_a_study_with_a_kinship_is_being_written() -> None:
    """A mixed model that asks for the approximation is refused, and is not
    given the exact test in silence.

    The approximation is the plan `gwas-logistic`, being one item over both
    mixed models. A study that made the exact test of every variant and
    reported that it had approximated nothing would give the user no way to
    tell that what they asked for did not happen.
    """
    phenotypes = _phenotypes()

    with pytest.raises(ValueError, match="approximation is being written"):
        calc_gwas(
            open_vcf(PANEL),
            phenotypes["cont"],
            TraitType.CONTINUOUS,
            covariates=phenotypes[["cov1", "cov2"]],
            kinship=_the_kinship_of_the_panel(),
            use_grammar_gamma_approx=True,
        )


def test_the_wald_test_is_the_linear_models_own_and_the_score_test_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """`test="wald"` is the test a linear model makes, and asking for it
    gives the study that asking for nothing gives.

    What a user reads in `result.test` is what they may write back into the
    call: a `calc_gwas` that refused `test` refused the value it had just
    given them. The score test is the one a linear model has not, and the
    core is what says so.
    """
    asked = _the_worked_example(worked_example, test=popnei.TestType.WALD)
    by_default = _the_worked_example(worked_example)

    assert by_default.test == popnei.TestType.WALD
    assert asked.test == popnei.TestType.WALD
    numpy.testing.assert_allclose(
        asked.stats["beta"].to_numpy(),
        by_default.stats["beta"].to_numpy(),
        rtol=0,
        atol=0,
    )
    with pytest.raises(
        ValueError, match="only test is the t test of the effect it fitted"
    ):
        _the_worked_example(worked_example, test="score")


def test_a_test_of_another_name_is_refused_with_the_two_names(
    worked_example: pathlib.Path,
) -> None:
    """A name that is of neither test is refused with both of them, which
    the core holds beside the two names of a trait."""
    with pytest.raises(ValueError, match="`wald`.*`score`.*`rao`"):
        _the_worked_example(worked_example, test="rao")


def test_a_test_that_is_no_name_is_refused_by_its_type(
    worked_example: pathlib.Path,
) -> None:
    """What is written for `test` and is not a name is refused by its type,
    as a trait that is no name is."""
    with pytest.raises(TypeError, match="`test` is"):
        _the_worked_example(worked_example, test=1)


def test_a_covariate_named_twice_is_refused(
    worked_example: pathlib.Path,
) -> None:
    """Two covariates of one name are refused, naming the two columns.

    The effects come back under the names of the columns of the design, so
    two of one name would be one entry of `covariate_effects` and a user who
    asked for that name would read one of the two without knowing which. A
    frame takes a repeated column label, and the two columns here hold
    different values, so nothing else of the study refuses them: the design
    has three independent columns.
    """
    covariates = pandas.DataFrame(
        numpy.column_stack([WORKED_EXAMPLE_COVARIATE, [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]]),
        index=list(WORKED_EXAMPLE_INDIVIDUALS),
        columns=["cov", "cov"],
    )

    with pytest.raises(ValueError, match="the covariate 'cov' is named twice"):
        _the_worked_example(worked_example, covariates=covariates)


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


def test_a_covariate_that_is_not_finite_is_refused_by_its_name(
    worked_example: pathlib.Path,
) -> None:
    """An infinity is a number to pandas and not to a fit, and it is refused
    in this layer, which has the name of the covariate and of the individual.

    It is what a covariate that came out of the user's own arithmetic as an
    infinity gives, a division by 0 among the causes. The core refuses it
    too, by the column of the design and the place of the individual, which
    is what a caller of `popnei._core` reads; neither message names a file,
    since what is wrong is wrong whatever variants are read.
    """
    covariates = _the_covariate()
    covariates.loc["i1", "cov"] = numpy.inf

    with pytest.raises(ValueError, match="covariate 'cov' at 'i1' is inf") as refused:
        _the_worked_example(worked_example, covariates=covariates)
    assert str(refused.value).startswith("the value of the covariate")


def test_a_phenotype_that_is_not_finite_is_refused_by_its_individual(
    worked_example: pathlib.Path,
) -> None:
    """An infinite trait is refused where the individual has a name.

    A missing phenotype is an individual that is not tested, and an infinity
    is not that: it would carry through the null model into the effect of
    every variant, so the user is told which individual to leave out.
    """
    trait = _the_trait()
    trait["i4"] = numpy.inf

    with pytest.raises(ValueError, match="the phenotype of 'i4' is inf"):
        _the_worked_example(worked_example, phenotype=trait)


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


def test_a_frequency_of_half_called_genotypes_can_pass_a_half(
    tmp_path: pathlib.Path,
) -> None:
    """`allele_freq` is a mean dosage over the ploidy and not the frequency
    of the allele the dosages count, so it can pass a half.

    The case is the spec's, in "How it is verified" of "The linear model":
    five individuals at `0/. 0/. 0/. 0/. 1/1`. The major allele is the one
    that is most frequent among the called **alleles**, which counts the
    called half of a half called genotype, so it is `0`, on four halves
    against two; the mean that becomes `allele_freq` is over the whole
    called **genotypes**, and `1/1` is the only one, so the mean dosage is 2
    and the frequency is 1.0. Every dosage of the variant is that mean, so
    it has no variance and no answer either.

    It is a divergence from plink2 and not from pyNei, which agrees; neither
    reference panel shows it, one having every genotype called and the other
    having them missing whole.
    """
    path = tmp_path / "half_called.vcf"
    individuals = ("h0", "h1", "h2", "h3", "h4")
    path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.4",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                "\t".join(
                    [
                        "#CHROM",
                        "POS",
                        "ID",
                        "REF",
                        "ALT",
                        "QUAL",
                        "FILTER",
                        "INFO",
                        "FORMAT",
                        *individuals,
                    ]
                ),
                "chr1\t1000\thalf\tA\tT\t.\tPASS\t.\tGT\t0/.\t0/.\t0/.\t0/.\t1/1",
            ]
        )
        + "\n"
    )
    trait = pandas.Series([1.0, 2.0, 3.0, 4.0, 5.0], index=list(individuals))

    result = calc_gwas(open_vcf(path, only_passed=False), trait, "continuous")

    assert result.stats["allele_freq"][0] == 1.0
    assert numpy.isnan(result.stats["beta"][0])
    assert numpy.isnan(result.stats["se"][0])
    assert numpy.isnan(result.stats["p_value"][0])


def test_the_private_module_names_a_design_that_is_not_contiguous(
    worked_example: pathlib.Path,
) -> None:
    """What a user who calls `popnei._core` themselves reads.

    The package builds the design itself and it lies row after row, so
    nothing a user writes reaches this message. When it is read it names the
    argument and what makes an array the core can take, as it does for the
    table of a principal component analysis: the core reads the design one
    row of covariates for each individual, and an array that lies column
    after column holds the same numbers in another matrix.
    """
    variants = open_vcf(worked_example, only_passed=False)
    by_columns = numpy.asfortranarray(
        numpy.column_stack([numpy.ones(6), WORKED_EXAMPLE_COVARIATE])
    )

    with pytest.raises(ValueError, match="ascontiguousarray") as refusal:
        _core.calc_gwas(
            variants._source,
            numpy.arange(6, dtype=numpy.uint64),
            numpy.asarray(WORKED_EXAMPLE_TRAIT, dtype=numpy.float64),
            by_columns,
            "continuous",
            None,
            None,
            False,
            False,
            variants._steps,
        )

    assert "`design`" in str(refusal.value)


def _the_private_call(variants, **written) -> None:
    """`popnei._core.calc_gwas` over the worked example, with `written` over
    the arguments the package would have built.

    It is what a user who calls the private module themselves reaches, and
    the three kinship cases below are the only place the binding's own
    refusals of a kinship are read: the package builds the matrix itself,
    cut and contiguous, so nothing a user writes through `calc_gwas` gets
    there.
    """
    asked = {
        "individuals": numpy.arange(6, dtype=numpy.uint64),
        "phenotype": numpy.asarray(WORKED_EXAMPLE_TRAIT, dtype=numpy.float64),
        "design": numpy.ascontiguousarray(
            numpy.column_stack([numpy.ones(6), WORKED_EXAMPLE_COVARIATE])
        ),
        "kinship": None,
        "use_grammar_gamma_approx": False,
    }
    asked.update(written)
    _core.calc_gwas(
        variants._source,
        asked["individuals"],
        asked["phenotype"],
        asked["design"],
        "continuous",
        None,
        asked["kinship"],
        asked["use_grammar_gamma_approx"],
        False,
        variants._steps,
    )


def test_the_private_module_names_a_kinship_that_is_not_contiguous(
    worked_example: pathlib.Path,
) -> None:
    """A kinship that lies column after column is refused by the name of the
    argument, as the design is.

    `pyo3.md` asks for this test of every matrix that crosses, and the
    package cannot reach it: it indexes the frame's array, which gives one
    that lies row after row. A matrix that lies the other way round is its
    own transpose, which for a kinship is the same matrix, so what this
    refusal is really for is the strided view of another array that
    `as_slice` would take.
    """
    variants = open_vcf(worked_example, only_passed=False)
    by_columns = numpy.asfortranarray(numpy.eye(6))

    with pytest.raises(ValueError, match="ascontiguousarray") as refusal:
        _the_private_call(variants, kinship=by_columns)

    assert "`kinship`" in str(refusal.value)


def test_the_private_module_refuses_a_kinship_of_another_size(
    worked_example: pathlib.Path,
) -> None:
    """A kinship that is not one row and one column for each tested
    individual is a ``RuntimeError``, which is what the three buffers of a
    study that do not hold it are.

    The package builds the matrix from the individuals it tested, so a size
    that does not fit them is a defect of popnei and not something a user
    wrote. The message says both sizes.
    """
    variants = open_vcf(worked_example, only_passed=False)

    with pytest.raises(RuntimeError, match="the kinship holds 25 values"):
        _the_private_call(variants, kinship=numpy.ascontiguousarray(numpy.eye(5)))


def test_the_private_module_names_the_cell_of_a_kinship_that_is_not_finite(
    worked_example: pathlib.Path,
) -> None:
    """A kinship with a value that is not finite is refused by the cell it is
    in, before any model is fitted.

    `Kinship.__post_init__` refuses such a matrix, so what reaches here is a
    frame written into afterwards and a caller of the private module. Left
    in, it would spread through the eigendecomposition into every eigenvalue
    and the study would come back with a NaN for every variant.
    """
    variants = open_vcf(worked_example, only_passed=False)
    with_a_nan = numpy.ascontiguousarray(numpy.eye(6))
    with_a_nan[2, 3] = numpy.inf
    with_a_nan[3, 2] = numpy.inf

    with pytest.raises(ValueError, match="row 2 and the column 3 of the kinship"):
        _the_private_call(variants, kinship=with_a_nan)


# What GMMAT 1.5.0's `glmmkin` fitted for the panel with every genotype
# called, from `tests/reference/gwas/gmmat.null_models.tsv`, which the
# reference script writes at full precision: the variance of the random
# effect of the kinship, which GMMAT calls `tau`, what is left over, which it
# calls `sigma2`, and the effects of the intercept, of `cov1` and of `cov2`.
OF_GMMAT_NULL = {
    "genetic_variance": 1.22161667529699,
    "residual_variance": 0.342359482266917,
    "intercept": 4.67802051309181,
    "cov1": 0.473360959469751,
    "cov2": 1.11027907093373,
}

# How far each of those five may be from GMMAT's: 1e-5 absolute, the bound of
# "How it is verified" of "The linear mixed model" of the spec, which is how
# far two restricted maximum likelihood searches land apart. Nothing of it is
# spent on rounding, the file being at full precision.
#
# Measured over the five on 24 September 2026: the worst is the genetic
# variance, 1.221e-6 away on Accelerate and 1.218e-6 on faer, 12 per cent of
# what is allowed, and the largest of the three effects is the 6.93e-7 of
# `cov2` on both. Where that 1.2e-6 comes from is the kinship and not the
# search: GMMAT was given the six printed digits of plink2's matrix and this
# test reads the float64 beside them. The bound is the spec's and is not
# lowered to two or three times what was measured, as the bounds on popnei's
# own arithmetic are.
OF_GMMAT_NULL_MODEL = 1e-5

# How far a `-log10(p_value)` of the Wald test may be from rrBLUP's: 1e-4
# absolute, from the same section, which asks it of all 1200 variants. It is
# the distance between two fits and not a printed digit, the file being at
# full precision, and 1e-4 in `-log10(p)` is 2.3e-4 of the p-value itself.
#
# Measured over the 1200 on 24 September 2026: the worst is `var0572`,
# 1.997e-5 on Accelerate and 1.996e-5 on faer, 20 per cent of what is
# allowed.
OF_RRBLUP = 1e-4

# How far `1 / se**2` of the score test may be from GMMAT's `VAR`, as a share
# of it, and how far a p-value may be from GMMAT's in `log10`: 1e-5 and 1e-4,
# both from the spec and both over all 1200 variants of both panels.
#
# The two files are printed to six significant digits, which rounds a value
# by up to 5e-6 of itself, so a share of the first bound goes on GMMAT's
# printing. `log10` shrinks a relative difference, so the second has more.
#
# Measured over the 1200 of each panel on 24 September 2026, the same to
# three digits on both backends: the worst `1 / se**2` is 4.43e-6 of `VAR` at
# `var0955` of the panel with every genotype called and 5.42e-6 at `var1060`
# of the panel with genotypes missing, 54 per cent of what is allowed; the
# worst p-value is 4.62e-5 and 4.76e-5 in `log10`, both at `var0185`, 48 per
# cent of what is allowed.
#
# How much of each of those is popnei's, measured on 24 September 2026. At
# `var0955` GMMAT prints 11.6504, so half of its last digit is 4.2917e-6 of
# the value against the 4.4268e-6 measured: at most 1.35e-7 of that
# difference is popnei's arithmetic. The comparison is not for that reason
# empty. The largest printing contribution anywhere in the file is 4.9098e-6,
# so the 1e-5 bound still fails if popnei's variance is about 5.7e-6 relative
# off, which is forty times above the signal it can resolve; and 282 of the
# 1200 variants are further from GMMAT than their own printing accounts for,
# so there is real signal in the column. The p-value half is the tighter of
# the two: 4.62e-5 in `log10` against a printing floor of 2.2e-6 is
# twenty-one times above it.
OF_GMMAT_VARIANCE = 1e-5
OF_GMMAT_P_VALUE = 1e-4

# How far a number of a mixed model may be from pyNei's: 1.5e-7, for the two
# variances and the effects of the null model as a share of pyNei's, for a
# `beta` and an `se` as a share of the `se` of their variant, and for a
# p-value as the distance between the two in `log10`.
#
# It is the one comparison of the spec that does not take the common 1e-9
# relative, and "How it is verified" of "The linear mixed model" says why:
# the criterion of the restricted maximum likelihood search is flat at its
# minimum, so an eigenvalue moving in its last bits moves the ratio of the
# two variances by about the square root of that, and every number of this
# model is built from that ratio. The two backends alone put the genetic
# variance of this build 3.155e-9 apart on the same kinship, 2.583e-9 of it,
# so a bound at 1e-9 would sit below the noise of the search and would pass
# or fail by rounding. The 9.7e-9 the spec gives, and this comment gave
# until 24 September 2026, is the cargo build's figure, which is not this
# one: the two were measured on different builds and the wrong one was
# copied here.
#
# So it was lowered until it failed, on both backends, over all 1200 variants
# of both panels under both tests, on 24 September 2026. It breaks at 5.17e-8
# on Accelerate and 4.65e-8 on faer, at the residual variance of the null
# model, which is `delta` times the genetic one and so carries the whole
# error of the search; the worst `beta` is 4.21e-8 and 3.78e-8 of its `se`,
# at the panel with genotypes missing under the score test, and the worst
# p-value is 4.32e-8 and 3.88e-8 in `log10`. This is 2.9 times the worst of
# those, and 58 times the 2.583e-9 the two backends sit apart, so it is a
# bound the comparison can fail.
OF_PYNEI_WITH_A_KINSHIP = 1.5e-7

# How many of the 5 causal variants have to be among the 10 smallest p-values
# of the mixed model, which is what "How it is verified" of "What every model
# shares" asks: at least 3. Measured on 24 September 2026, 4 of the 5 are, on
# both backends.
CAUSAL_AMONG_THE_SMALLEST = 3
THE_SMALLEST_P_VALUES = 10


def _the_kinship_of_the_panel() -> Kinship:
    """The kinship that `plink2 --make-rel square bin` wrote for the panel
    with every genotype called, at full precision.

    It is the one the mixed models are given, for both panels, as "How it is
    verified" of "What every model shares" of the spec asks: a kinship that
    came from neither popnei nor pyNei, and the one
    `tests/reference/gwas/make_reference.py` gave GMMAT and rrBLUP. GMMAT
    fits one null model with it and scores the variants of each panel against
    that fit, so the panel with genotypes missing is tested with the other
    panel's kinship here as it is there.
    """
    with gzip.open(
        REFERENCE_KINSHIP_DIR / "panel_called.plink2.rel.bin.gz", "rb"
    ) as stored:
        entries = numpy.frombuffer(stored.read(), dtype="<f8")
    lines = (
        (REFERENCE_KINSHIP_DIR / "panel_called.plink2.rel.id").read_text().splitlines()
    )
    individuals = [line for line in lines[1:] if line]
    assert len(individuals) == PANEL_NUM_INDIVIDUALS
    matrix = entries.reshape(PANEL_NUM_INDIVIDUALS, PANEL_NUM_INDIVIDUALS).copy()
    return Kinship(
        matrix=pandas.DataFrame(matrix, index=individuals, columns=individuals),
        num_vars=PANEL_NUM_VARS,
    )


def _the_mixed_study_of(panel: pathlib.Path, test: str | None, with_cov1: bool = True):
    """The study of `panel` with the kinship plink2 wrote and the test
    `test`, `None` taking the default of the model.

    `with_cov1` says whether the continuous covariate goes in beside the
    intercept and the binary one. rrBLUP takes every fixed effect as a
    factor, so it was given `cov2` alone and popnei is run with that one
    covariate for its comparison; GMMAT was given both.
    """
    phenotypes = _phenotypes()
    columns = ["cov1", "cov2"] if with_cov1 else ["cov2"]
    return calc_gwas(
        open_vcf(panel),
        phenotypes["cont"],
        TraitType.CONTINUOUS,
        covariates=phenotypes[columns],
        kinship=_the_kinship_of_the_panel(),
        test=test,
    )


def test_the_null_of_the_panel_is_gmmats_variances_and_effects() -> None:
    """The five numbers of the null model against GMMAT's `glmmkin`, within
    1e-5 absolute.

    The search has to be reproduced step for step or the variances move, so a
    `genetic_variance` near but not at 1.221617 means the search and one far
    from it means the criterion or the clamp of the eigenvalues at 0. The
    cargo test of the same five asserts them at the fit; what this one adds
    is that they reach Python under the names `NullModel` gives them, where
    GMMAT calls the two variances `tau` and `sigma2`.

    The `heritability` is asserted with them because it is the only number of
    `NullModel` built from the two variances rather than read off the fit,
    and because it is `None` for every other model.
    """
    of_gmmat = pandas.read_csv(
        REFERENCE_GWAS_DIR / "gmmat.null_models.tsv", sep="\t"
    ).set_index("model")

    result = _the_mixed_study_of(PANEL, "score")

    null = result.null_model
    assert null.model == GWASModel.LMM
    assert result.test == "score"
    assert null.num_individuals == PANEL_NUM_INDIVIDUALS
    # The file is read as well, so that a number of it that moved away from
    # the literals above would fail here and not quietly widen the check.
    assert dict(
        zip(
            OF_GMMAT_NULL,
            of_gmmat.loc["lmm", ["tau", "sigma2", "intercept", "cov1", "cov2"]],
            strict=True,
        )
    ) == pytest.approx(OF_GMMAT_NULL, rel=0, abs=0)
    found = {
        "genetic_variance": null.genetic_variance,
        "residual_variance": null.residual_variance,
        "intercept": null.covariate_effects["intercept"],
        "cov1": null.covariate_effects["cov1"],
        "cov2": null.covariate_effects["cov2"],
    }
    for name, expected in OF_GMMAT_NULL.items():
        assert found[name] == pytest.approx(expected, rel=0, abs=OF_GMMAT_NULL_MODEL), (
            f"{name} is {found[name]} and GMMAT gives {expected}"
        )
    genetic = OF_GMMAT_NULL["genetic_variance"]
    residual = OF_GMMAT_NULL["residual_variance"]
    assert null.heritability == pytest.approx(
        genetic / (genetic + residual), rel=0, abs=OF_GMMAT_NULL_MODEL
    )


def test_every_variant_of_the_panel_is_rrblups_wald_test() -> None:
    """The `-log10(p_value)` of all 1200 variants against rrBLUP's `GWAS`
    with `P3D = TRUE`, within 1e-4 absolute.

    `-log10(p)` is all rrBLUP reports, so it is what is compared. The study
    is run with `cov2` alone: rrBLUP takes every fixed effect as a factor, so
    the reference script gave it that covariate and no other, and a run that
    gave it both would get numbers that are close and not equal, which reads
    like a tolerance that is too tight and is not.

    The markers are compared first, so that every row is matched to the
    variant rrBLUP wrote it for and not to the row at the same place.
    """
    of_rrblup = pandas.read_csv(
        REFERENCE_GWAS_DIR / "rrblup.panel_called.lmm.tsv", sep="\t"
    )

    result = _the_mixed_study_of(PANEL, "wald", with_cov1=False)

    assert result.null_model.model == GWASModel.LMM
    assert result.test == "wald"
    assert len(of_rrblup.index) == PANEL_NUM_VARS
    assert list(result.stats["id"]) == list(of_rrblup["marker"])
    # rrBLUP writes the column under the name of the trait it was given.
    theirs = of_rrblup["cont"].to_numpy()
    ours = -numpy.log10(result.stats["p_value"].to_numpy())
    difference = numpy.abs(ours - theirs)
    worst = int(numpy.argmax(difference))
    assert difference[worst] <= OF_RRBLUP, (
        f"-log10(p) of {of_rrblup['marker'][worst]} is {ours[worst]} and "
        f"rrBLUP gives {theirs[worst]}, {difference[worst]} away against the "
        f"{OF_RRBLUP} allowed"
    )


@pytest.mark.parametrize("panel", PANELS)
def test_every_variant_of_a_panel_is_gmmats_score_test(panel: pathlib.Path) -> None:
    """The variance of the score and the p-value of all 1200 variants of both
    panels against GMMAT's `glmm.score`.

    GMMAT reports the variance of the score, which is `x' p x`, the
    denominator both tests are built on and what popnei gives as
    `1 / se**2`, and the p-value. The p-values span 23 orders of magnitude
    and what a user reads is the exponent, so they are compared in `log10`.

    The second panel is what says that a missing genotype takes the mean
    dosage of its variant, which GMMAT calls `impute2mean`: 3 in 100 of its
    genotypes are missing whole, and it is scored against the null model
    fitted with the kinship of the panel where none is, as the reference
    script scores it.
    """
    name = "panel_called" if panel == PANEL else "panel"
    of_gmmat = pandas.read_csv(
        REFERENCE_GWAS_DIR / f"gmmat.{name}.lmm.score.tsv", sep="\t"
    )

    result = _the_mixed_study_of(panel, "score")

    assert len(of_gmmat.index) == PANEL_NUM_VARS
    assert list(result.stats["id"]) == list(of_gmmat["SNP"])
    se = result.stats["se"].to_numpy()
    ours = 1.0 / (se * se)
    theirs = of_gmmat["VAR"].to_numpy()
    share = numpy.abs(ours - theirs) / theirs
    worst = int(numpy.argmax(share))
    assert share[worst] <= OF_GMMAT_VARIANCE, (
        f"1 / se**2 of {of_gmmat['SNP'][worst]} of {name} is {ours[worst]} "
        f"and GMMAT gives {theirs[worst]}, which is {share[worst]} of it "
        f"against the {OF_GMMAT_VARIANCE} allowed"
    )
    ours_p = result.stats["p_value"].to_numpy()
    theirs_p = of_gmmat["PVAL"].to_numpy()
    apart = numpy.abs(numpy.log10(ours_p / theirs_p))
    worst = int(numpy.argmax(apart))
    assert apart[worst] <= OF_GMMAT_P_VALUE, (
        f"the p-value of {of_gmmat['SNP'][worst]} of {name} is "
        f"{ours_p[worst]} and GMMAT gives {theirs_p[worst]}, {apart[worst]} "
        f"apart in log10 against the {OF_GMMAT_P_VALUE} allowed"
    )


@pytest.mark.parametrize("test", ["wald", "score"])
@pytest.mark.parametrize("panel", PANELS)
def test_every_variant_of_a_panel_with_a_kinship_is_pyneis(
    panel: pathlib.Path, test: str
) -> None:
    """Both libraries on the same VCF with the same kinship, over all 1200
    variants of each panel under each test.

    Both are given the kinship plink2 wrote, so what is compared is the
    search and the two tests and not two kinships. The variants that have no
    answer are asserted to be the same ones, which on both panels is none of
    them.

    The bound is `OF_PYNEI_WITH_A_KINSHIP` and not the `OF_PYNEI` the linear
    model is held to, seven orders of magnitude apart: every number of this
    model comes out of a search whose criterion is flat at its minimum, and
    the comment on that constant has where each backend broke.
    """
    phenotypes = _phenotypes()
    theirs = pynei_gwas(
        vars_from_vcf(panel),
        phenotypes["cont"],
        "continuous",
        covariates=phenotypes[["cov1", "cov2"]],
        kinship=PyneiKinship(
            matrix=_the_kinship_of_the_panel().matrix, num_vars=PANEL_NUM_VARS
        ),
        test=test,
    )

    ours = _the_mixed_study_of(panel, test)

    assert ours.individuals == tuple(theirs.samples)
    assert ours.null_model.model == theirs.null_model.model
    assert ours.test == theirs.test
    for what, found, expected in (
        (
            "the genetic variance",
            ours.null_model.genetic_variance,
            theirs.null_model.genetic_variance,
        ),
        (
            "the residual variance",
            ours.null_model.residual_variance,
            theirs.null_model.residual_variance,
        ),
    ):
        assert found == pytest.approx(expected, rel=OF_PYNEI_WITH_A_KINSHIP), what
    numpy.testing.assert_allclose(
        ours.null_model.covariate_effects.to_numpy(),
        theirs.null_model.covariate_effects.to_numpy(),
        rtol=OF_PYNEI_WITH_A_KINSHIP,
        atol=0,
    )
    numpy.testing.assert_allclose(
        ours.stats["allele_freq"].to_numpy(),
        theirs.stats["allele_freq"].to_numpy(),
        rtol=0,
        atol=OF_PLINK2_FREQUENCY,
    )
    se = theirs.stats["se"].to_numpy()
    for column in ("beta", "se"):
        difference = numpy.abs(
            ours.stats[column].to_numpy() - theirs.stats[column].to_numpy()
        )
        worst = numpy.nanmax(difference / se)
        assert worst <= OF_PYNEI_WITH_A_KINSHIP, (
            f"the worst {column} is {worst} of the `se` of its variant "
            f"against the {OF_PYNEI_WITH_A_KINSHIP} allowed"
        )
    ours_p = ours.stats["p_value"].to_numpy()
    theirs_p = theirs.stats["p_value"].to_numpy()
    assert (numpy.isnan(ours_p) == numpy.isnan(theirs_p)).all()
    numpy.testing.assert_allclose(
        numpy.log10(ours_p),
        numpy.log10(theirs_p),
        rtol=0,
        atol=OF_PYNEI_WITH_A_KINSHIP,
    )


def test_the_mixed_model_finds_the_variants_that_were_planted() -> None:
    """At least 3 of the 5 causal variants are among the 10 smallest p-values
    of the mixed model.

    The trait was simulated from the genotypes with five causal variants of
    effect 0.6 and a heritability of 0.5, and the three subpopulations differ
    in their mean so that the structure of the panel confounds it. This is
    the one test of the file that asks whether the study answers the question
    it is for rather than whether it answers it as another program does.
    Measured on 24 September 2026, 4 of the 5 are there on both backends.
    """
    causal = set(pandas.read_csv(REFERENCE_GWAS_DIR / "causal_vars.csv")["id"])

    result = _the_mixed_study_of(PANEL, "wald")

    assert len(causal) == 5
    smallest = result.stats.nsmallest(THE_SMALLEST_P_VALUES, "p_value")
    found = causal & set(smallest["id"])
    assert len(found) >= CAUSAL_AMONG_THE_SMALLEST, (
        f"{sorted(found)} of the causal variants {sorted(causal)} are among "
        f"the {THE_SMALLEST_P_VALUES} smallest p-values, "
        f"{list(smallest['id'])}"
    )


def test_the_two_tests_of_the_mixed_model_differ_in_the_error_alone() -> None:
    """The Wald test and the score test give the same `beta` for every
    variant and a different `se`, and the default of the model is the Wald
    test.

    Both are `num / den` with the same numerator and the same denominator,
    and what they differ in is how uncertain that effect is and which
    distribution it is read against. A study that answered two different
    effects would have its defect in the projection matrix and not in either
    test, which is what makes this the cheapest check of the two together.
    """
    wald = _the_mixed_study_of(PANEL, "wald")
    score = _the_mixed_study_of(PANEL, "score")
    by_default = _the_mixed_study_of(PANEL, None)

    assert wald.test == "wald"
    assert score.test == "score"
    assert by_default.test == "wald"
    assert list(by_default.stats["p_value"]) == list(wald.stats["p_value"])
    assert list(wald.stats["beta"]) == list(score.stats["beta"])
    differ = numpy.asarray(wald.stats["se"]) != numpy.asarray(score.stats["se"])
    assert differ.all()


def test_a_tested_individual_the_kinship_has_not_is_refused_by_name() -> None:
    """A kinship of half the panel, with a phenotype of all of it.

    The random effect is the relatedness of every pair that is tested, so
    there is nothing to put in the row of an individual the matrix has not,
    and the refusal names the first of them. The other way round is no error:
    a kinship of the whole panel with a phenotype of some of it is cut to
    those, which the next test asserts.
    """
    kinship = _the_kinship_of_the_panel()
    of_the_first_hundred = kinship.filter_individuals(kinship.individuals[:100])

    with pytest.raises(ValueError, match="'s100' is tested and is not one of"):
        calc_gwas(
            open_vcf(PANEL),
            _phenotypes()["cont"],
            TraitType.CONTINUOUS,
            kinship=of_the_first_hundred,
        )


def test_a_kinship_of_more_individuals_than_are_tested_is_cut_to_them() -> None:
    """A kinship of the whole panel with a phenotype of 100 of them gives
    what a kinship already cut to those 100 gives, to the bit.

    The kinship a user holds is of their panel, and a phenotype that leaves
    individuals out does not make it another matrix, so it is cut here as
    `Kinship.filter_individuals` cuts it. The cutting is the only thing that
    can differ between the two calls, and a study that read the uncut matrix
    would be fitted on the relatedness of individuals it never tested.
    """
    kinship = _the_kinship_of_the_panel()
    tested = list(kinship.individuals[:100])
    phenotypes = _phenotypes().loc[tested]

    of_the_panel = calc_gwas(
        open_vcf(PANEL),
        phenotypes["cont"],
        TraitType.CONTINUOUS,
        covariates=phenotypes[["cov1", "cov2"]],
        kinship=kinship,
    )
    of_the_hundred = calc_gwas(
        open_vcf(PANEL),
        phenotypes["cont"],
        TraitType.CONTINUOUS,
        covariates=phenotypes[["cov1", "cov2"]],
        kinship=kinship.filter_individuals(tested),
    )

    assert of_the_panel.individuals == tuple(tested)
    assert of_the_panel.null_model.num_individuals == 100
    assert list(of_the_panel.stats["beta"]) == list(of_the_hundred.stats["beta"])
    assert list(of_the_panel.stats["se"]) == list(of_the_hundred.stats["se"])


def test_the_kinship_is_read_in_the_order_the_source_has_the_individuals() -> None:
    """A kinship whose rows were given in another order gives the study the
    source's order gives.

    The phenotype, the rows of the design and the dosages of a block are read
    together row by row, and the kinship is the relatedness of those rows, so
    a matrix left in the order the user built it in would put one
    individual's relatedness against another's genotypes. No message can
    catch this one, both matrices being kinships of the same 200 individuals:
    what says it is that the two studies are equal.
    """
    kinship = _the_kinship_of_the_panel()
    backwards = kinship.filter_individuals(list(reversed(kinship.individuals)))
    phenotypes = _phenotypes()

    in_the_sources_order = _the_mixed_study_of(PANEL, "score")
    written_backwards = calc_gwas(
        open_vcf(PANEL),
        phenotypes["cont"],
        TraitType.CONTINUOUS,
        covariates=phenotypes[["cov1", "cov2"]],
        kinship=backwards,
        test="score",
    )

    assert list(written_backwards.stats["beta"]) == list(
        in_the_sources_order.stats["beta"]
    )
    assert list(written_backwards.stats["se"]) == list(in_the_sources_order.stats["se"])


def test_a_kinship_that_does_not_tell_the_variances_apart_gives_none_of_them(
    worked_example: pathlib.Path,
) -> None:
    """A study over an identity kinship gives its rows and no variances.

    The identity is what a user passes to mean no relatedness, and with it
    the model is the ordinary linear one whatever the split between the
    genetic variance and the residual one: the restricted maximum likelihood
    has nothing to choose between them and its criterion is flat over the
    whole grid, so which point wins is rounding. What the study gave before
    this, on 24 September 2026, was a ``heritability`` of 6.83e-05, which
    reads as a small number and is an arbitrary one; perturbing such a
    kinship by 1e-15 gave 6.5e-5, 7.1e-5 and 0.967 over three seeds.

    The three fields are ``None`` and not 0 and not NaN, which is what the
    linear model gives for the two it has not, so a user meets one way of
    saying that a number is not there. The rows are still there and still
    the worked example's, because the test is scale free.

    It is the meanwhile of **Open 3** of `docs/specs/gwas.md`, and the test
    of the panel above is its other end: that fit is not flat and reports
    all three.
    """
    result = _the_worked_example(
        worked_example, kinship=_the_kinship_of_the_worked_example()
    )

    assert result.null_model.model == GWASModel.LMM
    assert result.null_model.genetic_variance is None
    assert result.null_model.residual_variance is None
    assert result.null_model.heritability is None
    assert result.null_model.covariate_effects.index.tolist() == ["intercept", "cov"]
    assert result.null_model.num_individuals == 6
    for at, (variant, _freq, beta, _se, p_value) in enumerate(WORKED_EXAMPLE_ROWS):
        if numpy.isnan(beta):
            assert numpy.isnan(result.stats["beta"][at]), variant
            continue
        assert result.stats["beta"][at] == pytest.approx(
            beta, rel=OF_THE_WORKED_EXAMPLE
        ), variant
        # The p-value is the worked example's as well: with an identity
        # kinship the mixed model's Wald test is the linear model's t test,
        # and both have the 3 degrees of freedom of 6 individuals less 2
        # coefficients and the variant.
        assert result.stats["p_value"][at] == pytest.approx(
            p_value, rel=OF_THE_WORKED_EXAMPLE
        ), variant


def test_the_panel_tells_the_two_variances_apart_and_gives_all_three() -> None:
    """The panel's own kinship is not flat, so its study reports the two
    variances and the heritability.

    It is the other end of the test above: a rule that gave ``None``
    wherever a mixed model was fitted would pass that one and take the
    heritability away from every user. The numbers themselves are GMMAT's
    and are asserted where the null model is.
    """
    result = _the_mixed_study_of(PANEL, "score")

    assert result.null_model.genetic_variance is not None
    assert result.null_model.residual_variance is not None
    assert result.null_model.heritability is not None


def test_a_kinship_that_is_not_a_kinship_is_refused_by_its_type() -> None:
    """The matrix itself is what a user gives instead, and what that gave
    before this check was the `AttributeError` of an object with no `matrix`
    in it."""
    with pytest.raises(
        TypeError, match="`kinship` is a DataFrame, and the kinship a mixed"
    ):
        calc_gwas(
            open_vcf(PANEL),
            _phenotypes()["cont"],
            TraitType.CONTINUOUS,
            kinship=_the_kinship_of_the_panel().matrix,
        )


def _the_kinship_of_the_worked_example(without: str | None = None) -> Kinship:
    """The kinship of the six individuals of the worked example, or of the
    five that are not `without`: the identity, which is the relatedness of
    individuals with no recent ancestor in common.

    The calls it is written for are refused before any model is fitted, so
    what the matrix holds only has to be a kinship.
    """
    names = [
        name
        for name in WORKED_EXAMPLE_INDIVIDUALS
        if without is None or name != without
    ]
    return Kinship(
        matrix=pandas.DataFrame(numpy.eye(len(names)), index=names, columns=names),
        num_vars=3,
    )


def _a_kinship_written_into() -> Kinship:
    """A kinship of the worked example with one cell of a pair written into
    after it was built, which is what the core is there to catch.

    ``Kinship.__post_init__`` refuses a matrix that is not symmetric, and the
    frame it holds is mutable, so the check it made says nothing about what a
    study is given later. The eigendecomposition reads the lower triangle
    alone, so such a matrix was being read as that half mirrored.
    """
    kinship = _the_kinship_of_the_worked_example()
    kinship.matrix.iloc[0, 1] = 0.5
    return kinship


def _the_calls_that_are_refused(
    worked_example: pathlib.Path,
) -> dict:
    """The call of each case of `refusals_of_both_layers.json`, over the
    worked example, as a function that makes it.

    The TypeScript suite holds the same cases under the same names, and each
    suite is what writes the call in its own language: the names of the
    arguments differ between the two.
    """
    covariates = _the_covariate()
    a_copy = covariates.assign(twice=covariates["cov"])
    not_finite = covariates.copy()
    not_finite.loc["i2", "cov"] = numpy.inf
    of_five = covariates.drop(index="i5")
    of_three = _the_trait()
    of_three[["i3", "i4", "i5"]] = numpy.nan
    named = _the_trait()
    named["i9"] = 5.0
    a_name = _the_trait().astype(object)
    a_name["i2"] = "tall"
    empty = _the_trait().astype(object)
    empty["i2"] = ""
    covariate_of_names = covariates.astype(object)
    covariate_of_names.loc["i2", "cov"] = "north"
    an_infinity = _the_trait()
    an_infinity["i2"] = numpy.inf
    a_missing_value = covariates.copy()
    a_missing_value.loc["i2", "cov"] = numpy.nan
    of_one_value = pandas.Series(
        4.0, index=list(WORKED_EXAMPLE_INDIVIDUALS), dtype=float
    )
    binomial = pandas.Series(
        [0.0, 1.0, 0.0, 1.0, 0.0, 1.0], index=list(WORKED_EXAMPLE_INDIVIDUALS)
    )
    return {
        "a kinship that is not a kinship": lambda: _the_worked_example(
            worked_example, kinship="a matrix"
        ),
        "the grammar gamma approximation with no kinship": lambda: _the_worked_example(
            worked_example, use_grammar_gamma_approx=True
        ),
        "the grammar gamma approximation with a kinship": lambda: _the_worked_example(
            worked_example,
            kinship=_the_kinship_of_the_worked_example(),
            use_grammar_gamma_approx=True,
        ),
        "a tested individual the kinship has not": lambda: _the_worked_example(
            worked_example,
            kinship=_the_kinship_of_the_worked_example(without="i2"),
        ),
        "the score test": lambda: _the_worked_example(worked_example, test="score"),
        "a test of another name": lambda: _the_worked_example(
            worked_example, test="rao"
        ),
        "a trait of another name": lambda: _the_worked_example(
            worked_example, trait="quantitative"
        ),
        "a binomial trait": lambda: _the_worked_example(
            worked_example, phenotype=binomial, trait="binomial"
        ),
        "a covariate named intercept": lambda: _the_worked_example(
            worked_example, covariates=covariates.rename(columns={"cov": "intercept"})
        ),
        "a covariate that has no value for a tested individual": lambda: (
            _the_worked_example(worked_example, covariates=of_five)
        ),
        "a covariate that is not finite": lambda: _the_worked_example(
            worked_example, covariates=not_finite
        ),
        "a covariate that is a copy of another": lambda: _the_worked_example(
            worked_example, covariates=a_copy
        ),
        "an individual of the phenotype that the variants have not": lambda: (
            _the_worked_example(worked_example, phenotype=named)
        ),
        "fewer individuals than the design has columns plus two": lambda: (
            _the_worked_example(worked_example, phenotype=of_three)
        ),
        "a trait that is the same in every individual": lambda: _the_worked_example(
            worked_example, phenotype=of_one_value
        ),
        "a phenotype that is a name": lambda: _the_worked_example(
            worked_example, phenotype=a_name
        ),
        "a phenotype that is the empty string": lambda: _the_worked_example(
            worked_example, phenotype=empty
        ),
        "a covariate that is a name": lambda: _the_worked_example(
            worked_example, covariates=covariate_of_names
        ),
        "a phenotype that is an infinity": lambda: _the_worked_example(
            worked_example, phenotype=an_infinity
        ),
        "a covariate that has no value at an individual": lambda: _the_worked_example(
            worked_example, covariates=a_missing_value
        ),
        "the covariates explaining the whole of the trait": lambda: _the_worked_example(
            worked_example,
            covariates=_the_trait().to_frame(name="itself"),
            kinship=_the_kinship_of_the_worked_example(),
        ),
        "a kinship that is not symmetric": lambda: _the_worked_example(
            worked_example, kinship=_a_kinship_written_into()
        ),
        "an approximation that is not a boolean": lambda: _the_worked_example(
            worked_example, use_grammar_gamma_approx="no"
        ),
    }


def _the_calls_that_are_coerced(worked_example: pathlib.Path) -> dict:
    """The call of each case of the `coercions` of
    `refusals_of_both_layers.json`: the worked example with one of its
    values written as something that is not a number, which `float` reads
    the number out of.

    The TypeScript suite holds the same cases under the same names and
    asserts the same null model of each of them.
    """
    trait_of_strings = pandas.Series(
        [str(value) for value in WORKED_EXAMPLE_TRAIT],
        index=list(WORKED_EXAMPLE_INDIVIDUALS),
    )
    covariate_of_strings = pandas.DataFrame(
        {"cov": [str(value) for value in WORKED_EXAMPLE_COVARIATE]},
        index=list(WORKED_EXAMPLE_INDIVIDUALS),
    )
    # The covariate of the worked example is 0 and 1, which is what a
    # boolean is read as: False, True, False, True, False, True.
    covariate_of_booleans = pandas.DataFrame(
        {"cov": [value == 1.0 for value in WORKED_EXAMPLE_COVARIATE]},
        index=list(WORKED_EXAMPLE_INDIVIDUALS),
    )
    # The digits of any script are the digits `float` reads, so a trait
    # written with the full width ones of a spreadsheet is the same trait.
    in_full_width = str.maketrans("0123456789", "０１２３４５６７８９")
    trait_in_full_width = pandas.Series(
        [str(value).translate(in_full_width) for value in WORKED_EXAMPLE_TRAIT],
        index=list(WORKED_EXAMPLE_INDIVIDUALS),
    )
    return {
        "a phenotype written in full width digits": lambda: _the_worked_example(
            worked_example, phenotype=trait_in_full_width
        ),
        "a phenotype written as strings": lambda: _the_worked_example(
            worked_example, phenotype=trait_of_strings
        ),
        "a covariate written as strings": lambda: _the_worked_example(
            worked_example, covariates=covariate_of_strings
        ),
        "a covariate written as booleans": lambda: _the_worked_example(
            worked_example, covariates=covariate_of_booleans
        ),
    }


def _the_calls_with_no_phenotype_for_one(worked_example: pathlib.Path) -> dict:
    """The call of each case of the `no_phenotype_in_both_layers` of that
    same file: a value that both layers read as an individual with no
    phenotype, which is left untested.

    A value means no phenotype exactly where `float` of it gives NaN, which
    is NaN itself and the string `nan`. The second is the one a user does not
    expect, and it is the one a table of traits written by a program that
    prints NaN as text arrives with.
    """
    of_nan = _the_trait()
    of_nan["i2"] = numpy.nan
    of_the_word = _the_trait().astype(object)
    of_the_word["i2"] = "nan"
    return {
        "a phenotype that is NaN": lambda: _the_worked_example(
            worked_example, phenotype=of_nan
        ),
        "a phenotype that is the string nan": lambda: _the_worked_example(
            worked_example, phenotype=of_the_word
        ),
    }


def _the_calls_that_typescript_alone_refuses(worked_example: pathlib.Path) -> dict:
    """The call of each case of the `refused_in_typescript_alone` of that
    same file: a phenotype that says here, and not in TypeScript, that `i2`
    has none.

    `None` is what this layer has for the `null` and the `undefined` of
    TypeScript, so the two cases are the same call here: `pandas.isna` reads
    it as a missing value, which is what `dropna` drops in pyNei, while
    `float(None)` raises and is what TypeScript is written to.
    """
    of_none = _the_trait().astype(object)
    of_none["i2"] = None
    return {
        "a phenotype that is null": lambda: _the_worked_example(
            worked_example, phenotype=of_none
        ),
        "a phenotype that is undefined": lambda: _the_worked_example(
            worked_example, phenotype=of_none
        ),
    }


def _the_calls_of(which: str, calls: dict) -> list[tuple[dict, object]]:
    """Each case of the list `which` of `refusals_of_both_layers.json` with
    the call this suite writes for it.

    A case of the file this suite has no call for, and a call this suite
    makes that the file does not list, fail here: that is what binds the
    file to the two suites.
    """
    listed = json.loads(
        (REFERENCE_GWAS_DIR / "refusals_of_both_layers.json").read_text()
    )[which]

    assert listed, f"`{which}` of the file lists nothing"
    assert sorted(calls) == sorted(case["case"] for case in listed), (
        f"the cases of `{which}` of refusals_of_both_layers.json and the "
        f"calls this suite writes are not the same: a case of the file is "
        f"answered for in both suites, which is what makes it bind"
    )
    return [(case, calls[case["case"]]) for case in listed]


def test_both_layers_refuse_the_same_calls(worked_example: pathlib.Path) -> None:
    """Every call of the `refusals` of
    `tests/reference/gwas/refusals_of_both_layers.json` is a `ValueError`
    whose message holds what that file says.

    The TypeScript suite walks the same file, so a refusal that one layer
    has and the other has not is a case in the file that one of the two
    suites cannot make a call for, and that suite fails. Nothing else says
    that the two layers refuse the same things: each of them checks itself
    against its own copy of the literals, which is how they came to refuse
    `test` differently for a day.

    A case that says `python_raises` is one where what the user wrote is of
    the wrong type altogether, which this layer has an exception of its own
    for and which JavaScript has not: it is the class that differs and not
    the message, and the message is what the file binds.
    """
    for case, call in _the_calls_of(
        "refusals", _the_calls_that_are_refused(worked_example)
    ):
        raises = TypeError if case.get("python_raises") == "TypeError" else ValueError
        match = case.get("match", case.get("match_in_python"))
        with pytest.raises(raises, match=match):
            call()


def test_both_layers_read_a_value_that_is_not_a_number_as_the_number_it_holds(
    worked_example: pathlib.Path,
) -> None:
    """Every call of the `coercions` of that file gives the null model of
    the worked example.

    A trait read from a file arrives as strings often enough that refusing
    it would be a divergence users feel, so `float` here and `Number` in
    TypeScript read the number out of the value, which is what pyNei does
    and what the spec settled on 23 September 2026. Each of these calls is
    the worked example with one of its values written another way, so each
    of them gives the numbers the worked example gives.
    """
    for case, call in _the_calls_of(
        "coercions", _the_calls_that_are_coerced(worked_example)
    ):
        result = call()

        assert result.null_model.num_individuals == 6, case["case"]
        for name in ("intercept", "cov"):
            assert result.null_model.covariate_effects[name] == pytest.approx(
                WORKED_EXAMPLE_NULL[name], rel=OF_THE_WORKED_EXAMPLE
            ), f"the effect of {name} with {case['case']}"
        assert result.null_model.residual_variance == pytest.approx(
            WORKED_EXAMPLE_NULL["residual_variance"], rel=OF_THE_WORKED_EXAMPLE
        ), f"the residual variance with {case['case']}"


def test_both_layers_leave_an_individual_with_no_phenotype_untested(
    worked_example: pathlib.Path,
) -> None:
    """Every call of the `no_phenotype_in_both_layers` of that file leaves
    `i2` untested, here and in TypeScript.

    A value means no phenotype exactly where `float` of it gives NaN, which
    is the rule the spec settled on 23 September 2026 by the oracle, and the
    string `nan` is the case nobody guesses: `float('nan')` is NaN, so it is
    an individual that is not tested and not a value that is refused. The
    TypeScript suite asserts the same five individuals of the same calls.
    """
    for case, call in _the_calls_of(
        "no_phenotype_in_both_layers",
        _the_calls_with_no_phenotype_for_one(worked_example),
    ):
        result = call()

        assert result.individuals == ("i0", "i1", "i3", "i4", "i5"), case["case"]
        assert result.null_model.num_individuals == 5, case["case"]


def test_what_typescript_alone_refuses_is_an_individual_with_no_phenotype_here(
    worked_example: pathlib.Path,
) -> None:
    """Every call of the `refused_in_typescript_alone` of that file leaves
    `i2` untested here, where TypeScript refuses it.

    Python says that an individual has no phenotype with a name the series
    has not and with the missing value of pandas, which `None` is one
    spelling of and which `dropna` drops in pyNei. TypeScript has no such
    value, and `float(None)` raises, so `null` and `undefined` are refused
    there, naming the individual. The TypeScript suite asserts the other
    half of each of these two.
    """
    for case, call in _the_calls_of(
        "refused_in_typescript_alone",
        _the_calls_that_typescript_alone_refuses(worked_example),
    ):
        result = call()

        assert result.individuals == ("i0", "i1", "i3", "i4", "i5"), case["case"]
        assert result.null_model.num_individuals == 5, case["case"]


def test_what_is_not_a_variants_is_refused_by_its_type() -> None:
    """The path of the VCF is what a user gives instead, and what that gave
    before this check was the `AttributeError` of an object with no source
    inside it."""
    with pytest.raises(TypeError, match="`calc_gwas` reads the variants"):
        calc_gwas(str(PANEL), _the_trait(), "continuous")
