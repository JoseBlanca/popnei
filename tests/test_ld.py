"""The matrix of the r² of every pair of variants, from Python.

r² is the square of the correlation between the dosages of two variants,
the dosage of an individual being how many of the alleles of its genotype
are not the major allele of the variant, and it says how much the genotype
of one variant tells about the genotype of the other.
:func:`popnei.calc_rogers_huff_r2_matrix` gives it for every pair of the
variants of one pass, and :class:`popnei.R2Matrix` is what it comes in.

`docs/specs/ld.md` has the calculation and every literal of this file. The
reference program is plink2 v2.0.0-a.7.7, whose matrix of
`tests/reference/ld/ld.vcf.gz`, 500 variants of 100 diploid individuals, is
stored in that directory and is what the cargo tests of the core assert
against; the five pairs of the table of "How it is verified" are asserted
here again, through the three layers, within 1e-12 relative.

The comparison the spec asks for at this function is with pyNei at commit
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
change on either side of it shows.
"""

import math
import re
from pathlib import Path

import numpy
import pytest
from popnei import R2Matrix, _core, calc_rogers_huff_r2_matrix, open_vcf
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
    """popnei and pyNei over the variants of `ld.vcf.gz` that are called at
    every individual, each library running its own missing data filter.

    With no missing genotype the two rules are one rule, so the two matrices
    are the same calculation and agree within 1e-12 relative. pyNei gives r,
    which is squared before the comparison.
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
    what they are asking for. The cap is what the user wrote, which is wrong
    whatever file is read, so the message names no file.
    """
    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(_the_ld_dataset(), max_num_vars=100)

    message = str(refusal.value)
    # The pass is stopped at the block that passes the cap, so the variants
    # it names are the first count above it, which is the whole file when the
    # reader gives it in one block, and the memory is the square of that
    # count, 8 bytes a pair.
    said = re.match(r"the pass gave (\d+) variants and `max_num_vars` is 100", message)
    assert said is not None, message
    stopped_at = int(said.group(1))
    assert stopped_at > 100
    assert str(stopped_at * stopped_at * 8) in message
    assert not message.startswith(str(REFERENCE_LD_DIR))


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


def test_a_source_with_no_variant_is_refused(write_vcf) -> None:
    """A VCF whose header names three individuals and that has no data line.

    A matrix of no pair is no answer, so a pass that gives no variant is a
    wrong input, as it is for the distances of `docs/specs/dists.md`, and
    the message says whether the source had no variant or the steps kept
    none.
    """
    path = write_vcf([])

    with pytest.raises(ValueError) as refusal:
        calc_rogers_huff_r2_matrix(open_vcf(path))

    assert str(refusal.value) == (
        f"{path}: the source has no variant, and a calculation needs 1 variant at least"
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
        f"{path}: the steps kept no variant of the 2 the source gave, and a "
        f"calculation needs 1 variant at least: the filter `missing_data` "
        f"was given 2 variants and kept 0"
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
