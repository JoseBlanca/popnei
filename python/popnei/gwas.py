"""The association study: which variants are associated with a trait.

A trait is one number per individual, a measurement such as the height of a
plant or the 0 and 1 of an individual that has a condition and one that has
not. :func:`calc_gwas` fits a model of the trait once with no variant in it,
which is the null model, and then tests every variant of a
:class:`popnei.Variants` against what that model left unexplained: the effect
of one more copy of a non major allele, how uncertain that effect is, and the
probability of seeing an effect that far from 0 when the variant has none.

Covariates are other numbers per individual whose effect on the trait is not
of interest but has to be taken out, the field a plant grew in for instance,
and the top principal components of the panel go in as covariates so that a
variant which only marks ancestry does not look associated.

Three of the four models are built. A trait that is a measurement without a
kinship is the linear model, which is what plink2's ``--glm`` computes, and
its test is the t test of the effect: the effect divided by its standard
error, which under the hypothesis that the variant has none follows a
Student t distribution with as many degrees of freedom as there are
individuals left once the covariates and the variant have been fitted, and
the p-value is the chance that such a t falls further from 0 than this one
did, either way.

With a kinship it is the linear mixed model, which a panel with families in
it needs: the trait carries a random effect whose covariance is the kinship
times a variance, so that two related individuals are expected to resemble
each other before any variant is looked at. Its two tests are rrBLUP's Wald
test, which estimates the scale of the two variances again with the variant
in, and GMMAT's score test, which holds both at the null.

A trait that is 0 and 1 without a kinship is the logistic regression, which
is what plink2's ``--glm`` computes for such a trait: the chance that an
individual is a 1 is a logistic curve in the covariates and the variant, and
the effect is a log odds ratio. Its default is the Wald test, one fit per
variant with the variant in it, and it also takes the score test, which fits
nothing per variant and which R's ``anova(glm, test = "Rao")`` computes. A
variant that separates the individuals that have the condition from those
that have not has no finite effect, and its Wald test gives NaN for all
three numbers.

The logistic mixed model, a trait that is 0 and 1 with a kinship, is being
written, and so is the GRAMMAR-Gamma approximation that a mixed model can
take instead of the exact denominator of its test; asking for either is a
``ValueError`` that says so.

`docs/specs/gwas.md` has the four models, the numbers the tests assert and
what popnei does differently from pyNei.
"""

import math
from dataclasses import dataclass
from enum import StrEnum

import numpy
import pandas

from popnei import _core
from popnei.kinship import Kinship
from popnei.variant import PassStats, Variants, _pass_stats_of

# The name the effect of the column of ones every design has comes back
# under in `NullModel.covariate_effects`, which is pyNei's name for the same
# row of its series. A covariate of that name is refused, since the two
# would be one row of it.
_INTERCEPT = "intercept"


class TraitType(StrEnum):
    """What was measured on each individual."""

    CONTINUOUS = "continuous"
    """A measurement, one number per individual, such as the height of a
    plant. Without a kinship it is fitted by a linear model and with one by a
    linear mixed model."""

    BINOMIAL = "binomial"
    """0 or 1: an individual that has a condition and one that has not.
    Without a kinship it is fitted by a logistic regression and with one by a
    logistic mixed model."""


class TestType(StrEnum):
    """Which test is made of every variant.

    Both ask whether the effect of the variant on the trait is 0, and under
    that they have the same distribution in large samples; they differ in
    what they cost.
    """

    WALD = "wald"
    """The model is fitted again with the variant in it, and the variant's
    effect is measured in its own standard errors away from 0. It costs a fit
    per variant, so it is the default where that fit is cheap."""

    SCORE = "score"
    """How steeply the fit would improve if the variant's effect were let off
    0, measured at the null model and against how uncertain that slope is. It
    costs no fit per variant."""


class GWASModel(StrEnum):
    """Which of the four models a study fits, which the trait and the
    kinship together decide."""

    LM = "lm"
    """A linear model: a continuous trait and no kinship. The trait is a
    straight line in the covariates and the variant."""

    LMM = "lmm"
    """A linear mixed model: a continuous trait with the kinship as the
    covariance of a random effect."""

    GLM = "glm"
    """A logistic regression: a binomial trait and no kinship."""

    GLMM = "glmm"
    """A logistic mixed model: a binomial trait with the kinship as the
    covariance of a random effect."""


@dataclass(frozen=True)
class NullModel:
    """The model a study fitted with no variant in it, which every variant
    was then tested against.

    It is pyNei's ``NullModel`` of ``pynei/gwas.py``, whose ``num_samples``
    is :attr:`num_individuals`, the name `docs/glossary.md` gives.
    """

    model: GWASModel
    """Which of the four models it is."""

    covariate_effects: pandas.Series
    """The effect of the column of ones, under the name ``intercept``, and of
    every covariate, under the name it was given, in the units of the trait
    for a continuous one and as a log odds ratio for a binomial one.

    A covariate named ``intercept`` is refused at the call, since the two
    would be one row of this series."""

    residual_variance: float | None
    """What the model left unexplained, and ``None`` for a binomial trait,
    whose variance is decided by its mean."""

    genetic_variance: float | None
    """The variance of the random effect of the kinship, and ``None`` without
    a kinship."""

    heritability: float | None
    """The genetic variance over the sum of the two, only for the linear
    mixed model."""

    num_individuals: int
    """How many individuals the study tested."""


@dataclass(frozen=True)
class GWASResult:
    """What an association study gives back.

    It is pyNei's ``GWASResult`` of ``pynei/gwas.py``, with these
    differences: ``samples`` is :attr:`individuals`, the name
    `docs/glossary.md` gives; it carries :attr:`pass_stats`, the counts of
    the pass the study made, which pyNei keeps in its ``Variants``; and
    ``chrom``, ``pos`` and ``id`` are columns of :attr:`stats` whenever the
    source has them, where pyNei leaves each out when its chunks have no such
    column.
    """

    stats: pandas.DataFrame
    """One row for each variant the study was given, in the order the
    variants came.

    The columns are ``chrom``, ``pos`` and ``id`` when the source carries
    them, and then ``allele_freq``, ``beta``, ``se`` and ``p_value``.
    ``pos`` is unsigned, as ``Block.pos`` and ``R2Matrix.poss`` are, where
    pyNei's is signed: taking one position from another wraps instead of
    going below 0, so a user who wants the distance between two variants
    casts the column first, ``stats['pos'].astype('int64')``.
    ``allele_freq`` is the frequency of the alleles that are not the major
    one over the tested individuals, ``beta`` the effect of one more copy of
    such an allele, in the units of the trait for a continuous one and as a
    log odds ratio for a binomial one, ``se`` the standard error of that
    effect and ``p_value`` the probability of an effect that far from 0 when
    the variant has none. A variant whose dosages are all the same among
    the tested individuals has no variance and cannot be tested: its row is
    here with its ``allele_freq``, and the other three are NaN. So is a
    variant whose logistic fit does not settle, under the Wald test of a
    binomial trait alone: one that separates the individuals that have the
    condition from those that have not, whose effect has no finite value to
    walk towards, and one that repeats a covariate, which separates nobody
    and leaves the fit a system with no one solution. The score test fits
    nothing for a variant and gives all three numbers for either of them."""

    null_model: NullModel
    """The model fitted with no variant in it."""

    trait: TraitType
    """What was measured."""

    test: TestType
    """Which test was made of every variant."""

    individuals: tuple[str, ...]
    """The names of the individuals that were tested, those that have a
    phenotype, in the order the source has them, which is the order their
    phenotype and their design were read in."""

    used_grammar_gamma_approx: bool
    """Whether the GRAMMAR-Gamma approximation was used, which only a mixed
    model can use and which is being written.

    It stands in for the denominator of a mixed model's test, which is a
    product with the covariance of the random effect and costs one such
    product for every variant: the approximation computes one factor from
    the first variants of the pass and reuses it, which makes the test
    linear in the individuals per variant instead of quadratic, at the cost
    of accuracy where a panel is strongly structured."""

    pass_stats: PassStats
    """The counts of the pass the study made: how many variants it was given,
    after the steps of the ``Variants``, and how many each filter of it was
    given and kept."""


def calc_gwas(
    variants: Variants,
    phenotype: pandas.Series,
    trait: TraitType | str,
    covariates: pandas.DataFrame | None = None,
    kinship: Kinship | None = None,
    test: TestType | str | None = None,
    # The default is written here and not taken from the core's
    # `DEFAULT_USE_GRAMMAR_GAMMA_APPROX`, which is what a study that can
    # make the approximation takes when the user says nothing: popnei
    # refuses the approximation until it is written, so a core whose default
    # became true would turn every plain call into a refusal. Until then
    # this default changes no result, since the approximation is refused
    # whatever is written here.
    use_grammar_gamma_approx: bool = False,
    transform_to_biallelic: bool = _core.DEFAULT_TRANSFORM_TO_BIALLELIC,
) -> GWASResult:
    """Which of the variants of `variants` are associated with `phenotype`.

    Each variant becomes one number per individual, its dosage: how many
    alleles of the genotype are not the major allele of that variant, which
    is the most frequent among its called alleles and the lowest numbered of
    two that are equally frequent. A genotype with any allele missing takes
    the mean of the dosages of its variant. The dosages, that mean, the
    frequency and whether a variant varies at all are over the tested
    individuals and not over the whole panel, so a phenotype that leaves
    individuals out gives a variant another frequency than the panel's.

    `phenotype` is a series indexed by individual name, and a frame of one
    column is taken as that column. The individuals that are tested are those
    that have a value there which is not missing and that the `variants` has,
    in the order the `variants` has them, whatever order the phenotype was
    written in. A name of the phenotype that is of nobody the pass gives is a
    ``ValueError`` that names it, and so is a name that is there twice.

    A value of the phenotype that is not a number is read as one with
    ``float``, as pyNei reads it, so a trait written as the strings
    ``['2', '3', '5']`` and one written as ``True`` and ``False`` are both
    accepted. A value whose ``float`` is NaN, which the string ``'nan'`` is
    as much as NaN itself, is an individual with no phenotype and is left
    untested, and one whose ``float`` raises is a ``ValueError`` naming the
    individual. An infinity is refused too, which pyNei accepts and then
    gives NaN for every variant. A trait that is the same in every tested
    individual is a ``ValueError``: there is nothing for a variant to be
    associated with, and pyNei refuses such a trait only when it is
    binomial.

    `trait` is ``"continuous"``, a measurement, or ``"binomial"``, 0 for an
    individual that has not a condition and 1 for one that has, the two
    values of :class:`TraitType`. A binomial trait whose value at a tested
    individual is neither 0 nor 1 is a ``ValueError`` naming the place of
    that individual, and so is one where every tested individual has the
    same value, which leaves one of the two groups empty. Without a kinship
    a binomial trait is a logistic regression and ``beta`` is a log odds
    ratio; with one it is the logistic mixed model, which is being written,
    and asking for it is a ``ValueError`` that says so. A null model whose
    fit walks towards an infinite coefficient instead of settling is a
    ``ValueError`` too: what takes it there is a covariate that separates
    the individuals that have the condition from those that have not, and
    the user takes that covariate out.

    `covariates` is a frame indexed by individual with one column for each
    covariate, and the design of the study is a column of ones for the
    intercept and one column for each of them. A covariate that does not
    cover every tested individual, one that holds a value which is missing or
    is not a number, and one named ``intercept``, which is the name the
    effect of the column of ones comes back under, are each a ``ValueError``.
    A covariate whose values are names is given as one column for each of
    them, 1 for the individuals of that value and 0 for the others, which is
    what ``pandas.get_dummies`` writes. Covariates that are not independent,
    one that is constant or a copy of another, are a ``ValueError`` too, at
    the tolerance numpy's ``matrix_rank`` uses, so a design popnei refuses is
    a design pyNei refuses; and so is a study of no more individuals than the
    columns of its design plus one, which would leave nothing to measure the
    uncertainty of a variant's effect from.

    `kinship` is the relatedness of every pair, what :func:`popnei.calc_kinship`
    gives or a :class:`popnei.Kinship` built from a matrix another program
    wrote, and it makes the study a linear mixed model: the trait carries a
    random effect of that covariance, so that a variant which only marks the
    ancestry of a panel does not look associated. It has to hold every
    individual that is tested, and a tested individual it has not is a
    ``ValueError`` naming them; the rows and the columns of the ones it holds
    over are left out, as :meth:`popnei.Kinship.filter_individuals` leaves
    them out. Without a kinship the structure of a panel is accounted for
    with the top principal components of
    :meth:`popnei.Kinship.principal_components` or of
    :func:`popnei.do_pca_from_variants` as covariates, which is enough for
    individuals that are not close relatives, and both at once is the Q+K
    model of a strongly subdivided panel.

    `test` is ``"wald"`` or ``"score"``, the two values of
    :class:`TestType`, and ``None`` takes the default of the model, which is
    the Wald test for the three models that are built. The linear model's
    only test is the t test of the effect it fitted, which is the Wald test,
    so ``"score"`` is a ``ValueError`` that says so; the linear mixed model
    takes either, the Wald test being rrBLUP's and the score test GMMAT's,
    and so does the logistic regression, whose Wald test fits one logistic
    regression per variant and whose score test fits none.

    `use_grammar_gamma_approx` stands in for the denominator of a mixed
    model's test, which costs a product with the covariance of the random
    effect for every variant. It is being written, and asking for it is a
    ``ValueError``: with no kinship because there is no such denominator to
    approximate, and with one because popnei cannot approximate it yet.

    `transform_to_biallelic` makes every allele that is not the major one
    count the same, which is what a variant of more than two different
    alleles among its called genotypes needs: without it such a variant is a
    ``ValueError`` that says which one it is, as it is in
    :func:`popnei.calc_kinship` and :func:`popnei.do_pca_from_variants`.

    The call makes one pass over the source of `variants`, through the steps
    that are on it, so the study is over the variants its filters kept, and
    the ``Variants`` is as it was afterwards. A pass that gives no variant is
    a ``ValueError`` whose message says whether the source held none or the
    steps kept none.

    What it gives is a :class:`GWASResult` with one row per variant in
    ``stats``, the null model, the individuals that were tested and the
    counts of the pass.

    It is pyNei's ``calc_gwas`` of ``pynei/gwas.py``, with these differences:
    `samples` is `individuals` in the result; there is no `num_threads`,
    because no calculation of popnei has one; `transform_to_biallelic` is
    new, where pyNei collapses a multiallelic variant silently; the result
    carries the counts of the pass; and ``chrom``, ``pos`` and ``id`` are
    columns of ``stats`` whenever the source has them.
    """
    if not isinstance(variants, Variants):
        # What a user gives instead is usually the path of the VCF, and what
        # that gave was the `AttributeError` of an object with no source
        # inside it.
        raise TypeError(
            f"`variants` is {variants!r}, a {type(variants).__name__}, and "
            f"`calc_gwas` reads the variants of a source: give it what "
            f"`open_vcf` or `open_vars` gives, "
            f"calc_gwas(open_vcf(vcf_path), phenotype, 'continuous')"
        )
    _refuse_a_kinship_that_is_not_one(kinship)
    trait_name = _a_name_written_in(
        "trait",
        trait,
        "it says what was measured: write trait='continuous' for a "
        "measurement of each individual and trait='binomial' for 0 and 1, "
        "which are the two members of `TraitType`",
    )
    test_name = (
        None
        if test is None
        else _a_name_written_in(
            "test",
            test,
            "it says which test is made of every variant: write test='wald', "
            "which is the one a linear model has, or test='score', which a "
            "mixed model has too, or leave it out for the default of the "
            "model, and the two of them are the members of `TestType`",
        )
    )
    tested = _the_tested_individuals(_the_phenotype(phenotype), variants.individuals)
    names = [name for name, _, _ in tested]
    covariate_names, columns = _the_covariates(covariates, names)
    null, stats, used_grammar_gamma_approx, counts = _core.calc_gwas(
        variants._source,
        numpy.fromiter(
            (position for _, position, _ in tested),
            dtype=numpy.uint64,
            count=len(tested),
        ),
        numpy.fromiter(
            (value for _, _, value in tested), dtype=numpy.float64, count=len(tested)
        ),
        _the_design(columns, len(names)),
        trait_name,
        test_name,
        _the_kinship_of_the_tested(kinship, names),
        use_grammar_gamma_approx,
        transform_to_biallelic,
        variants._steps,
    )
    model, test_name, effects, residual, genetic, heritability, num_individuals = null
    return GWASResult(
        stats=_the_stats(stats),
        null_model=NullModel(
            model=GWASModel(model),
            covariate_effects=_the_effects(effects, covariate_names),
            residual_variance=residual,
            genetic_variance=genetic,
            heritability=heritability,
            num_individuals=num_individuals,
        ),
        trait=TraitType(trait_name),
        test=TestType(test_name),
        individuals=tuple(names),
        used_grammar_gamma_approx=used_grammar_gamma_approx,
        pass_stats=_pass_stats_of(counts),
    )


def _refuse_a_kinship_that_is_not_one(kinship: Kinship | None) -> None:
    """What was written in `kinship`, refused unless it is a ``Kinship``.

    # Raises

    ``TypeError`` saying what a kinship is. What a user gives instead is
    usually the matrix itself, a frame or an array, and what that gave was
    the ``AttributeError`` of an object with no ``matrix`` in it.
    """
    if kinship is None or isinstance(kinship, Kinship):
        return
    # The value is named by its type and not written out: what a user gives
    # here is usually the matrix itself, and the repr of a frame of 200
    # individuals is the message and the frame.
    raise TypeError(
        f"`kinship` is a {type(kinship).__name__}, and the kinship a mixed "
        f"model takes is a `Kinship`: give it what `calc_kinship` gives, or "
        f"build one over the matrix another program wrote, "
        f"Kinship(matrix=frame, num_vars=num_vars)"
    )


def _the_kinship_of_the_tested(
    kinship: Kinship | None, tested: list[str]
) -> numpy.ndarray | None:
    """The kinship cut to the individuals of `tested`, in their order, as the
    core reads it: one row and one column for each of them, row after row.

    The core is given numbers and holds no name to cut a matrix by, so the
    cutting is here, where the names are. An individual the kinship holds
    over is left out, as `Kinship.filter_individuals` leaves it out: the
    kinship of a panel is the one a user has, and a phenotype that leaves
    individuals out does not make it another matrix.

    The rows and the columns are taken in one indexing of the array the
    frame holds, and not with `.loc` and then `to_numpy` and then
    `ascontiguousarray`: those are three copies of the matrix where the core
    needs one, because a frame of one dtype lies column after column and
    comes back from `to_numpy` the wrong way round for the core, which reads
    it row after row. Measured at 3000 individuals, a matrix of 72 MB, on 24
    September 2026: 144.4 MB at the peak and 72.4 MB held with the three, and
    72.2 MB at the peak and 72.0 MB held with this, which at the 10000
    individuals of
    `docs/objectives.md` is 800 MB of transient saved. Indexing with two
    arrays of positions gives a new array that already lies row after row,
    which is the layout the binding crate refuses anything else in.

    Whether the entries of the matrix are finite numbers and whether it is
    symmetric are not checked here: `Kinship.__post_init__` refuses a matrix
    that is neither, and the core refuses both again by the cell they are
    in, which is what catches a frame written into after it was built.

    # Raises

    ``ValueError`` when a tested individual is not one of the kinship's: the
    random effect of a mixed model is the relatedness of every pair that is
    tested, and there is nothing to put in the row of an individual the
    matrix has not.
    """
    if kinship is None:
        return None
    of_the_matrix = {name: row for row, name in enumerate(kinship.individuals)}
    rows = []
    for name in tested:
        row = of_the_matrix.get(name)
        if row is None:
            raise ValueError(
                f"{name!r} is tested and is not one of the "
                f"{len(kinship.individuals)} individuals of the `kinship`, "
                f"which holds the relatedness of every pair that is tested: "
                f"give a kinship of them, `calc_kinship(variants)` for "
                f"instance, or leave that individual out of the phenotype"
            )
        rows.append(row)
    of_the_pairs = numpy.asarray(rows, dtype=numpy.intp)
    return kinship.matrix.to_numpy(dtype=numpy.float64)[
        numpy.ix_(of_the_pairs, of_the_pairs)
    ]


def _a_name_written_in(argument: str, value: object, says: str) -> str:
    """The name a user wrote in `argument`, as the core reads it.

    Which names there are is the core's, which is where a name that is of
    none of them is refused, so what this makes sure of is that a name and
    not something else was written: `trait` has no default, and what a user
    writes there instead is usually the covariates.

    # Raises

    ``TypeError`` when what was written is not a string, with `says` telling
    what the argument is for. A member of a ``StrEnum`` is a string.
    """
    if not isinstance(value, str):
        raise TypeError(
            f"`{argument}` is {value!r}, a {type(value).__name__}, and {says}"
        )
    return str(value)


def _the_phenotype(phenotype: pandas.Series) -> pandas.Series:
    """The trait of the individuals as a series indexed by their names.

    A frame of one column is taken as that column, which is what a user gets
    from `frame[["cont"]]` and from a file read with one column kept.

    # Raises

    ``TypeError`` when `phenotype` is neither, and ``ValueError`` for a frame
    of any other number of columns, which names no one column to test.
    """
    if isinstance(phenotype, pandas.DataFrame):
        _, num_columns = phenotype.shape
        if num_columns != 1:
            raise ValueError(
                f"`phenotype` is a frame of {num_columns} columns, and a "
                f"study is of one trait: give the column of the trait, "
                f"frame['cont'], and the others as `covariates`"
            )
        return phenotype.iloc[:, 0]
    if not isinstance(phenotype, pandas.Series):
        raise TypeError(
            f"`phenotype` is {phenotype!r}, a {type(phenotype).__name__}, and "
            f"a trait is a pandas series of one number for each individual, "
            f"indexed by their names"
        )
    return phenotype


def _the_tested_individuals(
    phenotype: pandas.Series, of_the_pass: tuple[str, ...]
) -> list[tuple[str, int, float]]:
    """The individuals that are tested: those of `of_the_pass` that have a
    phenotype, each with where it is in the pass and what was measured on it.

    The order is the source's and not the phenotype's, because the phenotype,
    the rows of the design and the dosages of a variant are read together row
    by row: a study that took them in the order the phenotype was written in
    would measure one individual's trait against another's genotypes.

    # Raises

    ``ValueError`` when a name of the phenotype is of nobody the pass gives,
    when a name is there twice, and when a value of it holds no number. A
    value that means no phenotype is an individual that is not tested and no
    error, and what means that is what ``float`` reads as NaN: NaN itself,
    the missing value of pandas, and the string ``'nan'``, which is the one
    a user does not expect.
    """
    _refuse_an_individual_that_is_there_twice(list(phenotype.index))
    of_the_source = set(of_the_pass)
    for name in phenotype.index:
        if name not in of_the_source:
            raise ValueError(
                f"{name!r} has a phenotype and is not one of the "
                f"{len(of_the_pass)} individuals these variants give: the "
                f"individuals that are tested are those of the `Variants` "
                f"that have one"
            )
    measured = dict(zip(phenotype.index, phenotype, strict=True))
    tested: list[tuple[str, int, float]] = []
    for position, name in enumerate(of_the_pass):
        if name not in measured:
            continue
        value = measured[name]
        # A missing phenotype is an individual that is not tested, which is
        # what a name that is not in the series is: `None` and `pandas.NA`
        # are what a table read from a file holds where a value is blank.
        if pandas.isna(value):
            continue
        try:
            number = float(value)
        except TypeError, ValueError:
            raise ValueError(
                f"the phenotype of {name!r} is {value!r}, a "
                f"{type(value).__name__}, and a trait is a number: a trait "
                f"whose values are names is not a trait of a study, and one "
                f"that is 0 and 1 is written as those numbers with "
                f"trait='binomial'"
            ) from None
        # A value whose `float` is NaN is an individual with no phenotype,
        # which is what `pandas.isna` answered above for NaN itself and for
        # the missing value of a table; the string 'nan' arrives here
        # instead, and `float` makes the same NaN of it. The spec of the
        # study settled it on 23 September 2026, by the oracle: a value
        # means no phenotype exactly where `float` of it gives NaN.
        if math.isnan(number):
            continue
        # An infinity is refused here, where the individual has a name: it
        # would carry through the null model into the effect of every
        # variant, and the core, which refuses it too, has the place of the
        # individual and not its name.
        if not math.isfinite(number):
            raise ValueError(
                f"the phenotype of {name!r} is {number}, and a study is "
                f"fitted on numbers: leave that individual out of the "
                f"phenotype, which is what leaves it untested"
            )
        tested.append((name, position, number))
    return tested


def _the_covariates(
    covariates: pandas.DataFrame | None, tested: list[str]
) -> tuple[list[str], list[list[float]]]:
    """The names of the covariates and the value of each of them for each of
    the individuals of `tested`, in their order.

    # Raises

    ``TypeError`` when `covariates` is not a frame. ``ValueError`` when one
    of them is named ``intercept``, when the frame names an individual twice,
    when a covariate has no value for a tested individual, and when a value
    is missing or is not a number.
    """
    if covariates is None:
        return [], []
    if not isinstance(covariates, pandas.DataFrame):
        raise TypeError(
            f"`covariates` is {covariates!r}, a {type(covariates).__name__}, "
            f"and the covariates of a study are a pandas frame with one "
            f"column for each of them, indexed by the names of the "
            f"individuals"
        )
    names = [str(name) for name in covariates.columns]
    for name in names:
        if name == _INTERCEPT:
            raise ValueError(
                f"a covariate is named {_INTERCEPT!r}, which is the name the "
                f"effect of the column of ones every design has comes back "
                f"under in `null_model.covariate_effects`: the two would be "
                f"one row of that series, so give the covariate another name"
            )
    _refuse_a_covariate_that_is_there_twice(names)
    _refuse_an_individual_that_is_there_twice(list(covariates.index))
    of_the_frame = dict(
        zip(covariates.index, range(len(covariates.index)), strict=True)
    )
    values = covariates.to_numpy(dtype=object)
    columns = []
    for column, name in enumerate(names):
        columns.append(
            [
                _the_value_of_the_covariate(
                    name, individual, values, of_the_frame, column
                )
                for individual in tested
            ]
        )
    return names, columns


def _the_value_of_the_covariate(
    name: str,
    individual: str,
    values: numpy.ndarray,
    of_the_frame: dict,
    column: int,
) -> float:
    """The value of the covariate `name` at `individual`, as a number.

    # Raises

    ``ValueError`` when the individual is not in the frame, and when the
    value is missing or is not a number. A covariate whose values are names
    is refused here, and what a user does with one is to give one covariate
    for each of its values, which is what ``pandas.get_dummies`` writes.
    """
    row = of_the_frame.get(individual)
    if row is None:
        raise ValueError(
            f"the covariate {name!r} has no value for {individual!r}, and a "
            f"covariate holds one for every individual that is tested"
        )
    value = values[row, column]
    if pandas.isna(value):
        raise ValueError(
            f"the value of the covariate {name!r} at {individual!r} is "
            f"missing, and a covariate holds a number for every individual "
            f"that is tested: leave that individual out of the phenotype, or "
            f"fill the value in"
        )
    try:
        number = float(value)
    except TypeError, ValueError:
        raise ValueError(
            f"the value of the covariate {name!r} at {individual!r} is "
            f"{value!r}, a {type(value).__name__}, and a covariate is a "
            f"number: one whose values are names is given as one covariate "
            f"for each of them, 1 for the individuals of that value and 0 "
            f"for the others, which is what `pandas.get_dummies` writes"
        ) from None
    # An infinity is a number to pandas and not to a fit, and this is the
    # layer that has the name of the covariate and of the individual: the
    # core refuses it as well, by their places among the columns and the
    # rows, which is what a caller of `popnei._core` reads.
    if not math.isfinite(number):
        raise ValueError(
            f"the value of the covariate {name!r} at {individual!r} is "
            f"{number}, and a study is fitted on numbers: it would carry "
            f"through the null model into the effect of every variant"
        )
    return number


def _refuse_a_covariate_that_is_there_twice(names: list[str]) -> None:
    """The first covariate whose name another column has too, refused.

    # Raises

    ``ValueError`` naming it and the two columns it is in. The effects of
    the null model come back under the names of the columns of the design,
    so two covariates of one name would be one entry of
    ``covariate_effects``, and a user who asked for that name would read one
    of the two without knowing which. It is what a covariate named
    ``intercept`` is refused for, one column over.
    """
    first_at: dict = {}
    for column, name in enumerate(names):
        if name in first_at:
            raise ValueError(
                f"the covariate {name!r} is named twice, at the columns "
                f"{first_at[name]} and {column}, and the effect of every "
                f"covariate comes back under its name in "
                f"`null_model.covariate_effects`: the two would be one row "
                f"of that series, so give one of them another name"
            )
        first_at[name] = column


def _the_design(columns: list[list[float]], num_individuals: int) -> numpy.ndarray:
    """The design of the study: one row for each tested individual, the 1 of
    the intercept and then the value of each covariate.

    It is built here and not in the core, which takes the matrix as it is
    read, row after row: the names of the covariates and their frame are this
    layer's, and the core is given numbers.
    """
    design = numpy.empty((num_individuals, len(columns) + 1), dtype=numpy.float64)
    design[:, 0] = 1.0
    for column, values in enumerate(columns, start=1):
        design[:, column] = values
    return design


def _the_stats(stats: tuple) -> pandas.DataFrame:
    """The rows of the result as a frame, one row per variant.

    `copy=False`: the four columns of the variants are the arrays the core
    filled, 32 bytes for each variant and 32 MB for a million, and nothing
    else holds them.
    """
    chroms, poss, ids, allele_freq, beta, se, p_value = stats
    columns: dict = {}
    if chroms is not None:
        columns["chrom"] = chroms
    if poss is not None:
        columns["pos"] = poss
    if ids is not None:
        columns["id"] = ids
    columns["allele_freq"] = allele_freq
    columns["beta"] = beta
    columns["se"] = se
    columns["p_value"] = p_value
    return pandas.DataFrame(columns, copy=False)


def _the_effects(effects: numpy.ndarray, covariates: list[str]) -> pandas.Series:
    """The effect of each column of the design under its name: the intercept
    first, which is the column of ones, and then the covariates in the order
    they were given.

    # Raises

    ``RuntimeError`` when the model answered another number of effects than
    the design has columns, which is a defect of popnei: a user reports it
    instead of looking for what they typed wrong.
    """
    names = [_INTERCEPT, *covariates]
    if len(effects) != len(names):
        raise RuntimeError(
            f"the null model gave {len(effects)} effects and its design has "
            f"{len(names)} columns, which is a defect of popnei; please "
            f"report it"
        )
    return pandas.Series(effects, index=names, copy=False)


def _refuse_an_individual_that_is_there_twice(names: list) -> None:
    """The first name that is in `names` twice, refused.

    # Raises

    ``ValueError`` naming it and the two places it is at. The phenotype and
    the covariates hold one value for each individual, and a name that is
    there twice would give two of them where one was asked for, and would
    weigh twice in the null model and in every variant.
    """
    first_at: dict = {}
    for place, name in enumerate(names):
        if name in first_at:
            raise ValueError(
                f"the individual {name!r} is named twice, at the places "
                f"{first_at[name]} and {place}, and a study reads one "
                f"phenotype and one row of the design for each individual "
                f"that is tested"
            )
        first_at[name] = place
