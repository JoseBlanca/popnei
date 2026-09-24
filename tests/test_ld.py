"""The matrix of the r² of every pair of variants, and how r² falls off with
distance, from Python.

r² is the square of the correlation between the dosages of two variants,
the dosage of an individual being how many of the alleles of its genotype
are not the major allele of the variant, and it says how much the genotype
of one variant tells about the genotype of the other.
:func:`popnei.calc_rogers_huff_r2_matrix` gives it for every pair of the
variants of one pass, and :class:`popnei.R2Matrix` is what it comes in.
:func:`popnei.calc_ld_and_dist_per_pop` gives, for each population, the mean
r² of the pairs of each bin of the distance between their two variants, and
:class:`popnei.LdAndDistPerPop` is what that comes in.

`docs/specs/ld.md` has the calculation and every literal of this file. The
reference program is plink2 v2.0.0-a.7.7, whose matrix of
`tests/reference/ld/ld.vcf.gz`, 500 variants of 100 diploid individuals, is
stored in that directory and is what the cargo tests of the core assert
against; the five pairs of the table of "How it is verified" are asserted
here again, through the three layers, within 1e-12 relative.

The ten bins of the fall-off are the first table of that same part, of one
population of every individual, and their r² comes from plink2 too.

The comparison the spec asks for at the matrix is with pyNei at commit
ef0ca6e, which `pyproject.toml` names, and it has two halves, because the
two libraries do not read a missing genotype the same way: popnei leaves
the individual out of the pair its genotype is missing at, which is what
plink2 does, and pyNei leaves it in with a dosage of -1. With no missing
genotype left in the dataset the two agree within 1e-12 relative. With the
missing genotypes in they differ, and by how much is the table of "Missing
genotypes" of the spec, measured on `tests/reference/dists/panel.vcf.gz`,
200 individuals and 1200 variants with 3 in 100 genotypes missing: over the
1.4 million pairs off the diagonal of that panel, the r² of pyNei is
0.0037 from plink2's at the median, 0.047 at the 99th percentile and 0.194
at the largest, where popnei's is plink2's to the bit. So the difference
between popnei and pyNei on that panel is that table, and the test asserts
it instead of asserting that the two agree: the divergence is pinned, and a
change on either side of it shows. At the fall-off no value is compared with
pyNei, which gives a sample of pairs drawn with no seed; what is compared is
the set of variants each population keeps at its own major allele frequency.
"""

import math
import re
from pathlib import Path

import numpy
import pytest
from popnei import (
    R2Matrix,
    _core,
    calc_ld_and_dist_per_pop,
    calc_rogers_huff_r2_matrix,
    open_vcf,
)
from pynei import filter_by_maf as pynei_filter_by_maf
from pynei import filter_samples as pynei_filter_samples
from pynei import vars_from_vcf
from pynei.ld_calc import _calc_rogers_huff_r2
from pynei.var_filters import filter_by_missing_data as pynei_filter_by_missing_data

REFERENCE_LD_DIR = Path(__file__).parent / "reference" / "ld"
PANEL_VCF = Path(__file__).parent / "reference" / "dists" / "panel.vcf.gz"

# How many variants `tests/reference/ld/ld.vcf.gz` holds, how many of them
# have no variance, which is one dosage among their called genotypes, and
# how many are called at every individual, which is what
# `filter_by_missing_data(0)` keeps. The first two are "How it is verified"
# of `docs/specs/ld.md`; the third was counted with pyNei on 23 September
# 2026 and is asserted against pyNei's own filter below.
NUM_VARS_OF_THE_LD_DATASET = 500
VARS_WITH_NO_VARIANCE = 68
VARS_WITH_NO_MISSING_GENOTYPE = 21

# How many variants and individuals `tests/reference/dists/panel.vcf.gz`
# holds, from "How it is verified" of `docs/specs/ld.md`.
NUM_VARS_OF_THE_PANEL = 1200

# The five pairs of the table of "How it is verified" of `docs/specs/ld.md`,
# which plink2 v2.0.0-a.7.7 gave on 22 September 2026: the chromosome and
# the position of each of the two variants and their r². The last pair is on
# two chromosomes, which this function gives like any other.
THE_PAIRS_OF_THE_TABLE = [
    ("chr1", 1000, "chr1", 2000, 0.353466669239891),
    ("chr1", 1000, "chr1", 3000, 0.39849991080910563),
    ("chr1", 1000, "chr1", 11000, 0.24053784261608957),
    ("chr1", 1000, "chr1", 250000, 0.0256751927810228),
    ("chr1", 1000, "chr2", 1000, 0.008140034754693937),
]

# How close popnei and a stored number have to be, relative: the tolerance
# of "How it is verified" of `docs/specs/ld.md`. It is there for a version of
# plink2 that computes the expression in another order and not for popnei's
# own rounding, which on 22 September 2026 gave plink2's bits for every pair
# of the dataset, so a difference of 1e-13 here is something to look at.
TOLERANCE = 1e-12

# The three numbers of the table of "Missing genotypes" of
# `docs/specs/ld.md` for the rule of pyNei, a missing genotype left in the
# correlation as a dosage of -1: over the pairs off the diagonal of the
# panel its r² differs from plink2's, and so from popnei's, by this much at
# the median, at the 99th percentile and at the largest. Each is asserted to
# the digits the spec gives it in.
THE_DIVERGENCE_OF_PYNEI = ((0.0037, 1e-4), (0.047, 1e-3), (0.194, 1e-3))

# The bins of the three tables of "How it is verified" of
# `docs/specs/ld.md`: the distances from 1 to 250000 base pairs cut into ten
# of 25000, with the smallest and the largest distance of each, both
# included.
THE_BOUNDS_OF_THE_TEN_BINS = [
    (1, 25_000),
    (25_001, 50_000),
    (50_001, 75_000),
    (75_001, 100_000),
    (100_001, 125_000),
    (125_001, 150_000),
    (150_001, 175_000),
    (175_001, 200_000),
    (200_001, 225_000),
    (225_001, 250_000),
]

# The first table of "How it is verified" of `docs/specs/ld.md`, the one
# population of every one of the 100 individuals at a `max_allowed_maf` of
# 0.95: for each of the ten bins, how many pairs it holds, the mean of their
# r² and its standard deviation. Every value is the one that table prints,
# which `docs/reports/ld-method/bins.py` worked out from the r² plink2
# v2.0.0-a.7.7 gives for these individuals and these variants, and which
# `tests/reference/ld/ld.bins.txt` holds again.
THE_BINS_OF_EVERY_INDIVIDUAL = [
    (8744, 0.20767885551844031, 0.20568652974479465),
    (7815, 0.07890359176062511, 0.08625049214414023),
    (6846, 0.03508441302751711, 0.0411917416898793),
    (5962, 0.02056917580626069, 0.026923590343909974),
    (5140, 0.015026451851395499, 0.02060604430437979),
    (4168, 0.011542104404978385, 0.015425275975059542),
    (3308, 0.01145438215819949, 0.015306615529156098),
    (2447, 0.012095572873545887, 0.016610730726125223),
    (1481, 0.015365254218410632, 0.020746495170409326),
    (530, 0.013266303346602112, 0.017768417874071147),
]

# The two populations of the second and the third table of the same part:
# `pop_a` the individuals `i000` to `i049` of `ld.vcf.gz` and `pop_b` `i050`
# to `i099`.
THE_TWO_POPS = {
    "pop_a": [f"i{individual:03d}" for individual in range(50)],
    "pop_b": [f"i{individual:03d}" for individual in range(50, 100)],
}

# How many pairs each of the ten bins holds for each of the two populations
# and the mean of their r², from the second and the third table, which are
# of a `max_allowed_maf` of 0.8. Those two tables leave the standard
# deviations out to stay readable, and `tests/reference/ld/ld.bins.txt`
# holds them; the cargo tests of the core assert them.
THE_PAIRS_AND_MEANS_OF_THE_TWO_POPS = {
    "pop_a": [
        (7394, 0.22226316432228382),
        (6564, 0.09426862122352),
        (5648, 0.04565040582359873),
        (4918, 0.030489796800475328),
        (4304, 0.024750005446480792),
        (3540, 0.02220350204444288),
        (2872, 0.01770136991605654),
        (2137, 0.02022495044871507),
        (1240, 0.01792337843122636),
        (438, 0.020745833685396994),
    ],
    "pop_b": [
        (7625, 0.21935192592998345),
        (6779, 0.08778756463719926),
        (5968, 0.0442939962514918),
        (5189, 0.0321335996857634),
        (4473, 0.02676089802517462),
        (3567, 0.02136589829323258),
        (2823, 0.02578147369672746),
        (2086, 0.02294629377272581),
        (1275, 0.022448095913909734),
        (415, 0.016086351215632733),
    ],
}

# How many of the 500 variants of `ld.vcf.gz` each of the three tables
# keeps, from the same part: 432 at the `max_allowed_maf` of 0.95 of the
# first table, and 396 and 402 at the 0.8 of `pop_a` and of `pop_b`, worked
# out over the individuals of each population alone.
VARS_AT_THE_MAF_OF_THE_FIRST_TABLE = 432
VARS_OF_POP_A = 396
VARS_OF_POP_B = 402


def _the_ld_dataset():
    """The `Variants` of `tests/reference/ld/ld.vcf.gz`, read as plink2 read
    it: diploid, and with the variants that failed their FILTER among
    them."""
    return open_vcf(REFERENCE_LD_DIR / "ld.vcf.gz", only_passed=False)


def _the_variant_at(matrix: R2Matrix, chrom: str, pos: int) -> int:
    """Which row of the matrix is the variant at `pos` of `chrom`."""
    for row, (its_chrom, its_pos) in enumerate(
        zip(matrix.chroms, matrix.poss, strict=True)
    ):
        if its_chrom == chrom and its_pos == pos:
            return row
    raise AssertionError(f"the matrix has no variant at {chrom}:{pos}")


def _pynei_r2_of(path: Path, no_missing_genotype: bool = False) -> numpy.ndarray:
    """The r² of every pair of the variants of `path` as pyNei gives it.

    pyNei's `calc_rogers_huff_r2_matrix` gives r and not r², under
    "What pyNei does that is odd" of `docs/specs/ld.md`, so what comes back
    here is its value squared. Its guard on the major allele frequency is
    turned off: it refuses a whole dataset that holds a variant above 0.95,
    which popnei has `filter_by_maf` take out instead, and the 68 variants
    of the LD dataset that have no variance are all above it.
    """
    variants = vars_from_vcf(path)
    if no_missing_genotype:
        variants = pynei_filter_by_missing_data(variants, max_allowed_missing_rate=0)
    dosages = numpy.concatenate(
        [chunk.gts.to_012() for chunk in variants.iter_vars_chunks()], axis=0
    )
    return _calc_rogers_huff_r2(dosages, dosages, check_no_mafs_above=None) ** 2


def _pynei_vars_at_a_maf_of(individuals: list[str], max_allowed_maf: float) -> int:
    """How many variants of `ld.vcf.gz` pyNei keeps for `individuals` at
    that major allele frequency.

    pyNei has no `calc_ld_and_dist_per_pop` that counts them: it makes one
    pass per population, by putting `filter_samples` and `filter_by_maf`
    around the `Variants` it was given, so that is what is asked here.
    """
    variants = pynei_filter_samples(
        vars_from_vcf(REFERENCE_LD_DIR / "ld.vcf.gz"), individuals
    )
    kept = pynei_filter_by_maf(variants, max_allowed_maf)
    return sum(chunk.gts.gt_values.shape[0] for chunk in kept.iter_vars_chunks())


def _off_the_diagonal(matrix: numpy.ndarray) -> numpy.ndarray:
    """The cells of a square matrix that are not on its diagonal, which are
    the pairs of two different variants."""
    return matrix[~numpy.eye(matrix.shape[0], dtype=bool)]


def test_the_five_pairs_of_the_spec_are_the_ones_plink2_gives() -> None:
    """The matrix of `ld.vcf.gz` through the three layers.

    The five pairs are the table of "How it is verified" of
    `docs/specs/ld.md`, which the cargo tests of the core assert against the
    whole stored matrix of plink2; what this one adds is that the values
    reach Python at the row and the column of the right variant, since the
    chromosome and the position of each row are what a user finds a pair by.
    """
    matrix = calc_rogers_huff_r2_matrix(_the_ld_dataset())

    assert matrix.r2.shape == (NUM_VARS_OF_THE_LD_DATASET, NUM_VARS_OF_THE_LD_DATASET)
    assert matrix.r2.dtype == numpy.float64
    assert len(matrix.chroms) == NUM_VARS_OF_THE_LD_DATASET
    assert len(matrix.poss) == NUM_VARS_OF_THE_LD_DATASET
    assert set(matrix.chroms) == {"chr1", "chr2"}
    for chrom, pos, other_chrom, other_pos, r2 in THE_PAIRS_OF_THE_TABLE:
        row = _the_variant_at(matrix, chrom, pos)
        column = _the_variant_at(matrix, other_chrom, other_pos)
        assert math.isclose(matrix.r2[row, column], r2, rel_tol=TOLERANCE), (
            f"the pair {chrom}:{pos}, {other_chrom}:{other_pos}"
        )
        assert matrix.r2[column, row] == matrix.r2[row, column]


def test_a_variant_with_no_variance_has_nan_in_its_row_its_column_and_its_diagonal() -> (
    None
):
    """The 68 variants of `ld.vcf.gz` whose called genotypes hold one dosage.

    Such a variant says nothing about any other, itself included, so it has
    no r² anywhere and its diagonal cell is NaN and not 1. plink2 gives NaN
    there too.
    """
    matrix = calc_rogers_huff_r2_matrix(_the_ld_dataset())

    no_variance = numpy.isnan(numpy.diagonal(matrix.r2))
    assert int(no_variance.sum()) == VARS_WITH_NO_VARIANCE
    assert numpy.all(numpy.isnan(matrix.r2[no_variance, :]))
    assert numpy.all(numpy.isnan(matrix.r2[:, no_variance]))
    with_variance = ~no_variance
    assert numpy.all(numpy.diagonal(matrix.r2)[with_variance] == 1.0)


def test_the_matrix_and_the_positions_cannot_be_written_into() -> None:
    """A result holds the memory the core filled, and nothing writes into it
    afterwards: it is what the user asked for and not a buffer of theirs."""
    matrix = calc_rogers_huff_r2_matrix(_the_ld_dataset())

    assert not matrix.r2.flags.writeable
    assert not matrix.poss.flags.writeable
    with pytest.raises(ValueError):
        matrix.r2[0, 0] = 0.0


def test_the_counts_of_the_pass_are_in_the_result() -> None:
    """A pass with no filter takes every variant of the source, and the
    `Variants` is as it was afterwards: a second call gives the same
    matrix."""
    variants = _the_ld_dataset()

    matrix = calc_rogers_huff_r2_matrix(variants)

    assert matrix.pass_stats.num_vars == NUM_VARS_OF_THE_LD_DATASET
    assert matrix.pass_stats.filtering == {}
    assert variants.steps == ()
    again = calc_rogers_huff_r2_matrix(variants)
    assert numpy.array_equal(again.r2, matrix.r2, equal_nan=True)


def test_a_filter_on_the_variants_is_run_and_counted() -> None:
    """The filter that keeps the variants called at every individual.

    It is the one the comparison with pyNei is made over, and what it keeps
    of the 500 variants is 21: the matrix is of those alone, and the counts
    of the pass say what the filter was given and what it kept.
    """
    variants = _the_ld_dataset()
    variants.filter_by_missing_data(0)

    matrix = calc_rogers_huff_r2_matrix(variants)

    assert matrix.r2.shape == (
        VARS_WITH_NO_MISSING_GENOTYPE,
        VARS_WITH_NO_MISSING_GENOTYPE,
    )
    assert matrix.pass_stats.num_vars == VARS_WITH_NO_MISSING_GENOTYPE
    counts = matrix.pass_stats.filtering["missing_data"]
    assert counts.vars_processed == NUM_VARS_OF_THE_LD_DATASET
    assert counts.vars_kept == VARS_WITH_NO_MISSING_GENOTYPE


def test_the_matrix_is_the_r_squared_of_pynei_when_no_genotype_is_missing() -> None:
    """popnei and pyNei over the 21 variants of `ld.vcf.gz` that are called
    at every individual, each library running its own missing data filter.

    With no missing genotype the two rules for a missing genotype are one
    rule, so the two matrices are the same calculation, and what is left
    between them is the order the two sum in: popnei takes the six sums of a
    pair out of matrix products on the BLAS the build links, and pyNei
    centres the dosages and sums them with numpy. pyNei gives r, which is
    squared before the comparison.

    The two are 1.19e-13 apart, relative, at the pair they differ most at,
    against the 1e-12 they are asserted within: eight times of room.
    Measured on 23 September 2026 over the 256 pairs of those 21 variants
    that have an r² in both libraries, on the owner's Apple M5 Pro with the
    BLAS the native build links there, Accelerate. It does not move between
    runs, and the run with one rayon thread gives the same bits. So the
    tolerance is there for a machine whose BLAS sums in another order, or a
    release of pyNei that does, and not for popnei's own rounding, which is
    what "How it is verified" of `docs/specs/ld.md` says of its own 1e-12.
    Whoever meets this test red reads the largest difference off it: one
    near 1.19e-13 and above 1e-12 is the ground moving under the test, and
    one far above is popnei.
    """
    variants = _the_ld_dataset()
    variants.filter_by_missing_data(0)

    ours = calc_rogers_huff_r2_matrix(variants)
    theirs = _pynei_r2_of(REFERENCE_LD_DIR / "ld.vcf.gz", no_missing_genotype=True)

    assert theirs.shape == ours.r2.shape
    # The variants with no variance among the 21 have NaN in both libraries,
    # and a matrix of nothing else would pass this test whatever the numbers
    # were.
    assert numpy.isfinite(ours.r2).any()
    numpy.testing.assert_allclose(ours.r2, theirs, rtol=TOLERANCE, atol=0)


def test_the_matrix_differs_from_pyneis_by_the_table_of_the_spec() -> None:
    """popnei and pyNei over the panel of 200 individuals and 1200 variants
    with 3 in 100 genotypes missing, with those genotypes left in.

    The two libraries read a missing genotype differently, so this asserts
    the difference and not equality. popnei leaves the individual out of the
    pair, which is plink2's rule and gives plink2's bits on this panel, and
    pyNei leaves it in as a dosage of -1, which moves its r² from plink2's
    by the three numbers of the table of "Missing genotypes" of
    `docs/specs/ld.md`. Those three numbers are therefore the difference
    between the two libraries here, and a change on either side shows as one
    of them moving.

    The cap is the one call of this file that raises it above the 1200
    variants of the panel without taking the default, so that a cap a user
    wrote is known to be read and not only the one the core holds.
    """
    matrix = calc_rogers_huff_r2_matrix(open_vcf(PANEL_VCF), max_num_vars=2000)
    theirs = _pynei_r2_of(PANEL_VCF)

    assert matrix.r2.shape == (NUM_VARS_OF_THE_PANEL, NUM_VARS_OF_THE_PANEL)
    difference = _off_the_diagonal(numpy.abs(matrix.r2 - theirs))
    assert not numpy.any(numpy.isnan(difference))
    (median, of_the_median), (percentile, of_it), (largest, of_the_largest) = (
        THE_DIVERGENCE_OF_PYNEI
    )
    assert numpy.median(difference) == pytest.approx(median, abs=of_the_median)
    assert numpy.percentile(difference, 99) == pytest.approx(percentile, abs=of_it)
    assert difference.max() == pytest.approx(largest, abs=of_the_largest)


def test_the_two_libraries_differ_on_the_ld_dataset_with_its_missing_genotypes() -> (
    None
):
    """The same two rules over `ld.vcf.gz`, whose 1502 genotypes of 50000 are
    missing.

    The spec measured by how much the two differ on the panel and not on
    this dataset, so what is asserted here is what the spec states of it:
    that they do differ, at pairs whose variants both have an r². It is the
    guard against the two rules quietly becoming one.
    """
    ours = calc_rogers_huff_r2_matrix(_the_ld_dataset())
    theirs = _pynei_r2_of(REFERENCE_LD_DIR / "ld.vcf.gz")

    both = ~numpy.isnan(ours.r2) & ~numpy.isnan(theirs)
    assert both.sum() > 0
    assert not numpy.allclose(ours.r2[both], theirs[both], rtol=TOLERANCE, atol=0)


def test_a_max_num_vars_below_the_variants_of_the_pass_is_refused() -> None:
    """The 500 variants of `ld.vcf.gz` with a cap of 100.

    The matrix holds one r² for each pair, so it grows with the square of
    the variants, and a pass of more than the cap is refused instead of the
    machine being asked for the memory. The message carries the variants the
    pass had given when it was stopped, the cap, and the bytes the matrix of
    those variants would have held, so that a user who raises the cap knows
    what they are asking for. It names the file the pass was reading, as the
    dataset a principal component analysis is too large for does: the cap
    and the variants that file holds decide together whether the cap is
    passed, so a user who runs over a directory of VCFs is told which one
    went over it. The cap that no matrix could be held under, the test
    below, names no file, because that one is wrong before any file is
    opened.
    """
    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(_the_ld_dataset(), max_num_vars=100)

    message = str(refusal.value)
    assert message.startswith(f"{REFERENCE_LD_DIR / 'ld.vcf.gz'}: "), message
    # The pass is stopped at the block that passes the cap, so the variants
    # it names are the first count above it, which is the whole file when the
    # reader gives it in one block, and the memory is the square of that
    # count, 8 bytes a pair.
    said = re.search(r"the pass gave (\d+) variants and `max_num_vars` is 100", message)
    assert said is not None, message
    stopped_at = int(said.group(1))
    assert stopped_at > 100
    assert str(stopped_at * stopped_at * 8) in message


def test_a_max_num_vars_whose_matrix_is_not_counted_is_refused() -> None:
    """A cap of 2^63 variants, whose matrix holds more values than this
    machine counts.

    The matrix holds the variants squared, which is counted in a number of
    64 bits natively and of 32 in WebAssembly, where 65536 variants already
    pass it. The cap is looked at before the pass, so a cap no matrix could
    be held under is refused at the call and not after the file has been
    read, and, being what the user wrote, it names no file.
    """
    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(_the_ld_dataset(), max_num_vars=2**63)

    message = str(refusal.value)
    assert message.startswith("a `max_num_vars` of 9223372036854775808"), message
    assert not message.startswith(str(REFERENCE_LD_DIR))


def test_a_max_num_vars_of_no_variants_is_refused_at_the_call() -> None:
    """A cap of 0 variants, which is a matrix of no pair.

    The core takes the 0 and stops the pass at its first variant, so what a
    user would read is the message of a dataset too large for the cap, "the
    pass gave 1 variants and `max_num_vars` is 0 ... raise `max_num_vars` or
    filter the variants", where what they asked for is a matrix of no
    variant. The 0 is refused at the call instead, before the source is
    read, under the name the user wrote it in and with no file named, as
    `calcRogersHuffR2Matrix` of the TypeScript package refuses a
    `maxNumVars` of 0.
    """
    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(_the_ld_dataset(), max_num_vars=0)

    assert str(refusal.value) == (
        "`max_num_vars` is 0, and it says how many of something there are: a "
        "whole number of 1 or more that this machine can count"
    )


def test_a_source_with_no_variant_is_refused(write_vcf) -> None:
    """A VCF whose header names three individuals and that has no data line.

    A matrix of no pair is no answer, so a pass that gives no variant is a
    wrong input, as it is for the distances of `docs/specs/dists.md`, and
    the message says whether the source had no variant or the steps kept
    none. It is the one case every calculation over a pass raises, and the
    core is what writes the sentence, so it is the sentence the distances
    give.
    """
    path = write_vcf([])

    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(open_vcf(path))

    assert str(refusal.value) == (
        f"{path}: the pass gave no variant and its source holds none: a "
        f"statistic of a pass is calculated over the variants it gives"
    )


def test_steps_that_kept_no_variant_are_refused_with_what_the_filter_counted(
    write_vcf,
) -> None:
    """Two variants, each with a missing genotype, and a filter that keeps
    the variants with none.

    The counts of a pass that could not be finished are otherwise lost, so
    the message carries them: the filter was given the two variants of the
    source and kept neither.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t./.",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
        ]
    )
    variants = open_vcf(path)
    variants.filter_by_missing_data(0)

    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(variants)

    assert str(refusal.value) == (
        f"{path}: the pass gave no variant: its source gave 2 and the steps "
        f"kept none of them, the `missing_data` filter was given 2 and kept "
        f"0; a statistic of a pass is calculated over the variants it gives"
    )


def test_a_variant_with_no_called_genotype_gives_a_row_of_nan(write_vcf) -> None:
    """Three variants, one of them called in nobody.

    A variant that nobody was called at is called together with no other
    variant at any individual, so every pair it is in has no r² and its row,
    its column and its diagonal cell are NaN. It is an answer and not an
    error: the other variants have their r² in the same matrix.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t./.\t./.",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
        ]
    )

    matrix = calc_rogers_huff_r2_matrix(open_vcf(path))

    assert numpy.all(numpy.isnan(matrix.r2[1, :]))
    assert numpy.all(numpy.isnan(matrix.r2[:, 1]))
    assert matrix.r2[0, 2] == pytest.approx(1.0, rel=TOLERANCE)
    assert matrix.poss.tolist() == [10, 20, 30]
    assert matrix.chroms == ("chr1", "chr1", "chr1")


def test_what_is_no_variants_and_no_cap_is_refused(write_vcf) -> None:
    """The two arguments of the function, each given what it cannot be.

    The path of the VCF in the place of the `Variants` is the mistake that
    is easiest to make, and a cap that is no number of variants says nothing
    a user can mean: 2.5 variants, `True`, which Python would pass on as the
    number 1, and a negative number are refused under the name the user
    wrote them in.
    """
    path = write_vcf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"])

    with pytest.raises(TypeError, match="variants"):
        calc_rogers_huff_r2_matrix(path)
    with pytest.raises(TypeError, match="max_num_vars"):
        calc_rogers_huff_r2_matrix(open_vcf(path), max_num_vars=2.5)
    with pytest.raises(TypeError, match="max_num_vars"):
        calc_rogers_huff_r2_matrix(open_vcf(path), max_num_vars=True)
    with pytest.raises(ValueError, match="max_num_vars"):
        calc_rogers_huff_r2_matrix(open_vcf(path), max_num_vars=-1)


def test_the_default_cap_is_the_one_of_the_core() -> None:
    """The number of variants a user gets when they name none, which the
    core holds and the package puts in its signature: nothing writes 5000
    twice."""
    assert _core.DEFAULT_MAX_NUM_VARS == 5000


def test_ld_and_dist_gives_the_ten_bins_of_the_spec_that_plink2_gives() -> None:
    """The fall-off of r² with distance over `ld.vcf.gz`, through the three
    layers.

    The ten rows are the first table of "How it is verified" of
    `docs/specs/ld.md`: one population of every one of the 100 individuals
    at a `max_allowed_maf` of 0.95, which 432 of the 500 variants pass, with
    the distances from 1 to 250000 base pairs cut into ten bins of 25000.
    Every r² of every bin is plink2's, and the binning is the arithmetic the
    spec defines; `tests/reference/ld/ld.bins.txt` holds the same numbers
    with the standard deviations.

    The counts of pairs are compared exactly and the means and the standard
    deviations within 1e-12 relative, since both sides add the same r² and
    only the order of the sum can differ.
    """
    of_the_pass = calc_ld_and_dist_per_pop(
        _the_ld_dataset(),
        min_dist=1,
        max_dist=250_000,
        num_bins=10,
        max_allowed_maf=0.95,
    )

    # With no `pops` there is one population of every individual, named as
    # pyNei names it.
    assert list(of_the_pass.per_pop) == ["pop"]
    assert of_the_pass.num_vars_per_pop == {"pop": VARS_AT_THE_MAF_OF_THE_FIRST_TABLE}
    # The pass counted every variant of the file: the major allele frequency
    # takes variants out of a population and not out of the pass.
    assert of_the_pass.pass_stats.num_vars == NUM_VARS_OF_THE_LD_DATASET

    frame = of_the_pass.per_pop["pop"]
    assert list(frame.index) == [smallest for smallest, _ in THE_BOUNDS_OF_THE_TEN_BINS]
    assert frame.index.name == "smallest_dist"
    assert list(frame["largest_dist"]) == [
        largest for _, largest in THE_BOUNDS_OF_THE_TEN_BINS
    ]
    assert list(frame["num_pairs"]) == [
        num_pairs for num_pairs, _, _ in THE_BINS_OF_EVERY_INDIVIDUAL
    ]
    for row, (_, mean_r2, sd_r2) in enumerate(THE_BINS_OF_EVERY_INDIVIDUAL):
        assert frame["mean_r2"].iloc[row] == pytest.approx(mean_r2, rel=TOLERANCE)
        assert frame["sd_r2"].iloc[row] == pytest.approx(sd_r2, rel=TOLERANCE)


def test_ld_and_dist_gives_every_count_and_every_distance_as_a_signed_number() -> None:
    """The counts and the distances of the bins are signed 64 bit integers,
    as pyNei's counts are, so that the difference of two of them is a
    negative number and not 1.8e19.

    Over the one population of `ld.vcf.gz` at the settings of the first
    table of "How it is verified" of `docs/specs/ld.md`, the second bin
    holds 929 pairs fewer than the first, which an unsigned subtraction
    gives as 18446744073709550687.
    """
    of_the_pass = calc_ld_and_dist_per_pop(
        _the_ld_dataset(),
        min_dist=1,
        max_dist=250_000,
        num_bins=10,
        max_allowed_maf=0.95,
    )

    frame = of_the_pass.per_pop["pop"]
    num_pairs = frame["num_pairs"]
    assert int(num_pairs.iloc[1] - num_pairs.iloc[0]) == -929
    assert num_pairs.dtype == numpy.int64
    assert frame["largest_dist"].dtype == numpy.int64
    assert frame.index.dtype == numpy.int64


def test_ld_and_dist_keeps_in_each_pop_the_variants_pynei_keeps_there() -> None:
    """The variants `pop_a` and `pop_b` keep at a major allele frequency of
    0.8, counted by both libraries.

    "How it is verified" of `docs/specs/ld.md` compares no value with pyNei
    here, since pyNei gives a sample of pairs drawn with no seed, of r and
    not r², with a missing genotype left in. What is compared is the set of
    variants each population keeps: popnei counts them in one pass, in
    `num_vars_per_pop`, and pyNei is asked for them with `filter_samples`
    around the individuals of the population and `filter_by_maf` over what
    is left. Both have to give the 396 of `pop_a` and the 402 of `pop_b` of
    the spec.

    At 0.95 the two populations pass the same 432 variants, so 0.8 is the
    threshold that fails when the major allele frequency is worked out over
    all the individuals instead of over those of the population.
    """
    of_the_pass = calc_ld_and_dist_per_pop(
        _the_ld_dataset(),
        pops=THE_TWO_POPS,
        min_dist=1,
        max_dist=250_000,
        num_bins=10,
        max_allowed_maf=0.8,
    )

    # The populations come back in the order of the `pops` dict.
    assert list(of_the_pass.per_pop) == ["pop_a", "pop_b"]
    assert of_the_pass.num_vars_per_pop == {
        "pop_a": VARS_OF_POP_A,
        "pop_b": VARS_OF_POP_B,
    }
    assert of_the_pass.num_vars_per_pop == {
        pop: _pynei_vars_at_a_maf_of(individuals, 0.8)
        for pop, individuals in THE_TWO_POPS.items()
    }
    # The two populations count their own pairs, which the counts of the
    # first bin show: the variants they keep are not the same variants.
    assert list(of_the_pass.per_pop["pop_a"]["num_pairs"]) == [
        num_pairs for num_pairs, _ in THE_PAIRS_AND_MEANS_OF_THE_TWO_POPS["pop_a"]
    ]
    assert list(of_the_pass.per_pop["pop_b"]["num_pairs"]) == [
        num_pairs for num_pairs, _ in THE_PAIRS_AND_MEANS_OF_THE_TWO_POPS["pop_b"]
    ]
    for pop, rows in THE_PAIRS_AND_MEANS_OF_THE_TWO_POPS.items():
        means = of_the_pass.per_pop[pop]["mean_r2"]
        for row, (_, mean_r2) in enumerate(rows):
            assert means.iloc[row] == pytest.approx(mean_r2, rel=TOLERANCE)


def test_ld_and_dist_leaves_every_bin_empty_when_no_pair_reaches_min_dist(
    write_vcf,
) -> None:
    """A dataset of one chromosome whose variants span less than `min_dist`,
    which "The cases" of `docs/specs/ld.md` says is no error.

    Every bin holds 0 pairs and NaN for its mean and its standard deviation,
    and the three variants still passed the major allele frequency, so
    `num_vars_per_pop` counts them.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0",
        ]
    )

    of_the_pass = calc_ld_and_dist_per_pop(
        open_vcf(path), min_dist=1000, max_dist=5000, num_bins=4
    )

    frame = of_the_pass.per_pop["pop"]
    assert of_the_pass.num_vars_per_pop == {"pop": 3}
    assert list(frame.index) == [1000, 2001, 3001, 4001]
    assert list(frame["largest_dist"]) == [2000, 3000, 4000, 5000]
    assert list(frame["num_pairs"]) == [0, 0, 0, 0]
    assert frame["mean_r2"].isna().all()
    assert frame["sd_r2"].isna().all()


def test_ld_and_dist_refuses_what_is_no_variants_no_distance_and_no_bins(
    write_vcf,
) -> None:
    """The six arguments of the function, each given what it cannot be.

    The path of the VCF in the place of the `Variants` is the mistake that
    is easiest to make. A distance that is negative cannot reach the core,
    whose distances are unsigned, so the binding is what refuses it, and it
    is refused under the name the user wrote it in; so is one that is no
    whole number at all, 2.5 and `True`, which Python would pass on as the
    number 1. The other four are the core's: a `min_dist` above `max_dist`,
    a `num_bins` of 0, a `max_allowed_maf` outside 0 to 1, and a population
    that names an individual the dataset has not.
    """
    path = write_vcf(["chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"])

    with pytest.raises(TypeError, match="variants"):
        calc_ld_and_dist_per_pop(path)
    with pytest.raises(ValueError) as of_the_smallest_dist:
        calc_ld_and_dist_per_pop(open_vcf(path), min_dist=-1)
    # The smallest the two distances take is 0, which counts the pairs of
    # two variants at one position, and not the 1 of the window of
    # `filter_by_ld`: the refusal states the limit the argument has.
    assert "`min_dist` is -1" in str(of_the_smallest_dist.value)
    assert "0 or more" in str(of_the_smallest_dist.value)
    with pytest.raises(ValueError) as of_the_largest_dist:
        calc_ld_and_dist_per_pop(open_vcf(path), max_dist=-250)
    assert "`max_dist` is -250" in str(of_the_largest_dist.value)
    assert "0 or more" in str(of_the_largest_dist.value)
    with pytest.raises(ValueError, match=r"`num_bins` is -3"):
        calc_ld_and_dist_per_pop(open_vcf(path), num_bins=-3)
    with pytest.raises(TypeError, match="min_dist"):
        calc_ld_and_dist_per_pop(open_vcf(path), min_dist=2.5)
    with pytest.raises(TypeError, match="max_dist"):
        calc_ld_and_dist_per_pop(open_vcf(path), max_dist=True)
    with pytest.raises(TypeError, match="max_allowed_maf"):
        calc_ld_and_dist_per_pop(open_vcf(path), max_allowed_maf="a half")
    with pytest.raises(ValueError, match="min_dist"):
        calc_ld_and_dist_per_pop(open_vcf(path), min_dist=5000, max_dist=4000)
    with pytest.raises(ValueError, match="num_bins"):
        calc_ld_and_dist_per_pop(open_vcf(path), num_bins=0)
    with pytest.raises(ValueError, match="max_allowed_maf"):
        calc_ld_and_dist_per_pop(open_vcf(path), max_allowed_maf=1.5)
    with pytest.raises(ValueError, match="ind9"):
        calc_ld_and_dist_per_pop(open_vcf(path), pops={"pop1": ["ind1", "ind9"]})


def test_the_defaults_of_ld_and_dist_are_the_ones_of_the_core() -> None:
    """The four numbers a user gets when they name none, which the core
    holds and the package puts in its signature: nothing writes 1000000
    twice."""
    assert _core.DEFAULT_MIN_DIST == 1
    assert _core.DEFAULT_MAX_DIST == 1_000_000
    assert _core.DEFAULT_NUM_DIST_BINS == 50
    assert _core.DEFAULT_MAX_ALLOWED_MAF == 0.95
