"""The Kosman distances between individuals, from Python.

`docs/specs/dists.md` has the distance, the `Distances` it comes in and the
numbers this file asserts. Every literal here is the spec's or is read from
`tests/reference/dists/`, which
`tests/reference/dists/make_reference.py` writes: the four VCFs, the
distance of every pair that `gd.kosman` of PopGenReport 3.1.3 gives for each
of them, and the distance vector of pyNei for the two diploid ones.

The comparison the spec asks for is with pyNei at commit ef0ca6e, which
`pyproject.toml` names, and it is exact: both libraries divide the same two
integers once, so the vectors have to agree bit for bit, with NaN in the
same places. pyNei takes diploid genotypes alone, so the tetraploid and the
haploid datasets are checked against R's numbers and against the literals of
the spec's table.

The datasets that no file holds are written here as small VCFs: the four
tests of pyNei's `test/test_dists.py` take their genotypes, and pyNei is
given the same genotypes through its `Variants.from_gt_array`, which builds
a `Variants` of pyNei from an array.
"""

import math
from collections.abc import Sequence
from pathlib import Path

import numpy
import pandas
import pytest
from popnei import (
    Distances,
    FilteringStats,
    PassStats,
    calc_pairwise_kosman_dists,
    open_vars,
    open_vcf,
    write_vars,
)
from popnei.variant import Variants
from pynei import Variants as PyneiVariants
from pynei import calc_pairwise_kosman_dists as pynei_kosman_dists
from pynei import vars_from_vcf

REFERENCE_DISTS_DIR = Path(__file__).parent / "reference" / "dists"

# The four datasets of "How it is verified" of the spec, each with the ploidy
# its VCF is read with, how many variants and how many individuals it holds.
PANEL = ("panel", 2, 1200, 200)
FOUR_ALLELES = ("four_alleles", 2, 300, 40)
TETRAPLOID = ("tetraploid", 4, 200, 12)
HAPLOID = ("haploid", 1, 200, 12)

# The `min_num_snps` of the spec for the panel, which leaves the pairs that
# were called together at fewer than 1125 of its 1200 variants without a
# distance.
PANEL_MIN_NUM_SNPS = 1125

# The distances of the first three pairs of the tetraploid and the haploid
# datasets, t00 with t01, t02 and t03 and h00 with h01, h02 and h03, from the
# table of "How it is verified" of the spec, which `gd.kosman` gave. The
# tolerance is the spec's, the digits of the shortest of its numbers.
TETRAPLOID_LITERALS = [0.375, 0.38797814207650272, 0.40163934426229508]
HAPLOID_LITERALS = [
    0.62222222222222223,
    0.66847826086956519,
    0.63128491620111726,
]
TOLERANCE = 1e-9

# The diploid worked example of the spec, 4 variants of 3 individuals, as the
# data lines of a VCF and as the array of genotypes pyNei is given. The third
# variant has the half called genotype `0/.`, which is a missing genotype,
# and the fourth has a missing one.
WORKED_EXAMPLE_LINES = [
    "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t20\t.\tA\tT,C\t.\tPASS\t.\tGT\t0/1\t0/1\t1/2",
    "chr1\t30\t.\tA\tT,C\t.\tPASS\t.\tGT\t0/0\t0/.\t2/2",
    "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t1/1\t1/1",
]
WORKED_EXAMPLE_GTS = [
    [[0, 0], [0, 1], [1, 1]],
    [[0, 1], [0, 1], [1, 2]],
    [[0, 0], [0, -1], [2, 2]],
    [[-1, -1], [1, 1], [1, 1]],
]
# Its distances, the ploidy times the sum of d over the ploidy times n of the
# spec's table: 1 over 2 x 2, 5 over 2 x 3 and 2 over 2 x 3.
WORKED_EXAMPLE_DISTS = [1 / 4, 5 / 6, 1 / 3]
# The variants of each of its three pairs, which are the `min_num_snps` at
# which a pair keeps its distance and the one below which it loses it.
WORKED_EXAMPLE_NUM_VARS = [2, 3, 3]

# The genotypes of pyNei's `test_kosman_2_indis`, four individuals over 11
# variants, and the three distances it asserts: 1/3 for a with b, 0 for c
# with d, which hold the same genotype at every variant, and 0.45 for b with
# d. A genotype that is missing is missing whole, as it is in that test.
KOSMAN_2_INDIS_A = [
    [-1, -1],
    [0, 0],
    [0, 1],
    [0, 0],
    [0, 0],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 0],
    [0, 0],
    [0, 1],
]
KOSMAN_2_INDIS_B = [
    [1, 1],
    [-1, -1],
    [0, 0],
    [0, 0],
    [1, 1],
    [0, 1],
    [1, 0],
    [1, 0],
    [1, 0],
    [0, 1],
    [1, 1],
]
KOSMAN_2_INDIS_C = [[1, 1]] * 11
KOSMAN_2_INDIS_D = [[1, 1]] * 11
DIST_A_B = 1 / 3
DIST_C_D = 0.0
DIST_B_D = 0.45

# The genotypes of pyNei's `test_kosman_missing`: the same two individuals
# with one variant missing in each, and then with both of those variants
# missing in both of them. The two give the same distance, because a variant
# that either of them is not called at counts for neither.
MISSING_IN_EACH = (KOSMAN_2_INDIS_A, KOSMAN_2_INDIS_B)
MISSING_IN_BOTH = (
    [[-1, -1], [-1, -1]] + KOSMAN_2_INDIS_A[2:],
    [[-1, -1], [-1, -1]] + KOSMAN_2_INDIS_B[2:],
)

# The genotypes of pyNei's `test_kosman_pairwise`, which are those of
# `test_kosman_2_indis` with the last genotype of b holding the alleles 1 and
# 2, and the vector of six distances it asserts.
KOSMAN_PAIRWISE_B = KOSMAN_2_INDIS_B[:-1] + [[1, 2]]
KOSMAN_PAIRWISE_DISTS = [1 / 3, 0.75, 0.75, 0.5, 0.5, 0.0]

# The dataset of pyNei's `test_kosman_pairwise_with_filtered_vars`: 60
# variants of 30 individuals, every genotype called, which is what makes a
# missing data filter at 1 take nothing out.
FILTERED_VARS_NUM_VARS = 60
FILTERED_VARS_NUM_INDIVIDUALS = 30
FILTERED_VARS_SEED = 7

# The four individuals of "What pyNei does that is odd", with a distance that
# says which pair it is, and the lower triangle of their square matrix with
# its diagonal, which is what `triang_list_of_lists` gives.
ODD_VECTOR = [1, 2, 3, 12, 13, 23]
ODD_TRIANGLE = [[0], [1, 0], [2, 12, 0], [3, 13, 23, 0]]

# A vector of four distances, which is of no number of individuals: three
# individuals make three pairs and four make six.
VECTOR_OF_NO_NUMBER_OF_INDIVIDUALS = [0.1, 0.2, 0.3, 0.4]

_NUCLEOTIDES = "ACGT"


def _vcf_text(gts: numpy.ndarray, individuals: Sequence[str]) -> str:
    """The genotypes `gts`, variants x individuals x ploidy, as a VCF.

    Every variant declares the alleles of the whole dataset, one nucleotide
    each, and an allele below 0 is written as a dot: which letter an allele
    gets changes no distance, since the Kosman distance counts the alleles
    that two genotypes hold in common.
    """
    num_vars = gts.shape[0]
    alleles = _NUCLEOTIDES[: max(2, int(gts.max()) + 1)]
    lines = [
        "##fileformat=VCFv4.4",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + list(individuals)
        ),
    ]
    for index in range(num_vars):
        genotypes = [
            "/".join("." if allele < 0 else str(int(allele)) for allele in genotype)
            for genotype in gts[index]
        ]
        lines.append(
            "\t".join(
                [
                    "chr1",
                    str((index + 1) * 100),
                    ".",
                    alleles[0],
                    ",".join(alleles[1:]),
                    ".",
                    "PASS",
                    ".",
                    "GT",
                ]
                + genotypes
            )
        )
    return "\n".join(lines) + "\n"


@pytest.fixture
def vcf_of(tmp_path: Path):
    """A function that writes an array of genotypes as a VCF and gives its
    path, for the datasets that no file of `tests/reference/` holds."""

    def write(
        gts, individuals: Sequence[str] | None = None, name: str = "written.vcf"
    ) -> Path:
        gts = numpy.asarray(gts)
        if individuals is None:
            individuals = [f"i{index}" for index in range(gts.shape[1])]
        path = tmp_path / name
        path.write_text(_vcf_text(gts, individuals))
        return path

    return write


def _of_individuals(*individuals) -> numpy.ndarray:
    """The genotypes of a dataset whose individuals are given one by one,
    each of them the genotype it holds at every variant."""
    return numpy.transpose(numpy.array(individuals), axes=(1, 0, 2))


def _reference(name: str) -> Path:
    """The gzipped VCF of one of the four datasets of the spec."""
    return REFERENCE_DISTS_DIR / f"{name}.vcf.gz"


def _stored_dists(name: str, program: str) -> numpy.ndarray:
    """The distance of every pair that a program gave for a dataset, read
    from what the reference script stored: the first column of the file of
    `gd.kosman` and the only one of the file of pyNei."""
    numbers = numpy.loadtxt(
        REFERENCE_DISTS_DIR / f"{name}.{program}.tsv", skiprows=1, ndmin=2
    )
    return numbers[:, 0]


@pytest.mark.parametrize("dataset", [PANEL, FOUR_ALLELES])
def test_the_distances_of_a_diploid_dataset_are_the_ones_pynei_gives(
    dataset,
) -> None:
    """The two diploid datasets, read by both libraries.

    The vectors have to be equal bit for bit and not close: both libraries
    divide the same two whole numbers once, and a difference in the last
    place would say that one of them adds its distances in another way.
    """
    name, ploidy, _, num_individuals = dataset

    ours = calc_pairwise_kosman_dists(open_vcf(_reference(name), ploidy=ploidy))
    theirs = pynei_kosman_dists(vars_from_vcf(_reference(name)))

    assert ours.dist_vector.shape == (num_individuals * (num_individuals - 1) // 2,)
    assert numpy.array_equal(ours.dist_vector, theirs.dist_vector, equal_nan=True)
    assert numpy.array_equal(
        ours.dist_vector, _stored_dists(name, "pynei"), equal_nan=True
    )


def test_a_min_num_snps_leaves_the_pairs_pynei_leaves_without_a_distance() -> None:
    """The panel with `min_num_snps=1125`, in both libraries.

    Some pairs of the panel were called together at fewer than 1125 of its
    1200 variants, so the vector has NaN in it, and the two libraries have
    to have it in the same places.
    """
    name, ploidy, _, _ = PANEL

    ours = calc_pairwise_kosman_dists(
        open_vcf(_reference(name), ploidy=ploidy), min_num_snps=PANEL_MIN_NUM_SNPS
    )
    theirs = pynei_kosman_dists(
        vars_from_vcf(_reference(name)), min_num_snps=PANEL_MIN_NUM_SNPS
    )

    assert numpy.isnan(ours.dist_vector).any()
    assert numpy.array_equal(
        numpy.isnan(ours.dist_vector), numpy.isnan(theirs.dist_vector)
    )
    assert numpy.array_equal(ours.dist_vector, theirs.dist_vector, equal_nan=True)


@pytest.mark.parametrize(
    ("dataset", "literals"),
    [(TETRAPLOID, TETRAPLOID_LITERALS), (HAPLOID, HAPLOID_LITERALS)],
)
def test_the_distances_of_another_ploidy_are_the_ones_r_gives(
    dataset, literals
) -> None:
    """The tetraploid and the haploid datasets, which pyNei refuses and
    `gd.kosman` of R computes.

    The three literals of the spec's table are the first three pairs of each,
    and every pair of the dataset is compared with what R gave for it, so
    that the ploidy is read and not assumed to be 2 anywhere.
    """
    name, ploidy, _, num_individuals = dataset

    ours = calc_pairwise_kosman_dists(open_vcf(_reference(name), ploidy=ploidy))

    assert ours.dist_vector.shape == (num_individuals * (num_individuals - 1) // 2,)
    for index, expected in enumerate(literals):
        assert ours.dist_vector[index] == pytest.approx(expected, abs=TOLERANCE)
    assert numpy.allclose(
        ours.dist_vector, _stored_dists(name, "gdkosman"), atol=TOLERANCE, rtol=0
    )


def test_the_worked_example_gives_the_distances_of_the_spec_and_of_pynei(
    write_vcf,
) -> None:
    """The 4 variants of 3 individuals of "How it is verified", written as a
    VCF, with its half called genotype and its missing one.

    The three distances are the spec's integers divided as the spec divides
    them, and pyNei is given the same genotypes as an array.
    """
    ours = calc_pairwise_kosman_dists(open_vcf(write_vcf(WORKED_EXAMPLE_LINES)))
    theirs = pynei_kosman_dists(
        PyneiVariants.from_gt_array(
            numpy.array(WORKED_EXAMPLE_GTS), samples=["ind1", "ind2", "ind3"]
        )
    )

    assert list(ours.dist_vector) == WORKED_EXAMPLE_DISTS
    assert numpy.array_equal(ours.dist_vector, theirs.dist_vector, equal_nan=True)
    assert ours.names == ("ind1", "ind2", "ind3")


def test_a_pair_keeps_its_distance_at_exactly_min_num_snps_and_loses_it_below(
    write_vcf,
) -> None:
    """The worked example with `min_num_snps` of 3 and of 4.

    Its pairs were called together at 2, 3 and 3 variants, so 3 leaves the
    first pair without a distance and leaves the other two theirs, which is
    what "strictly below" means, and 4 leaves no pair a distance.
    """
    assert WORKED_EXAMPLE_NUM_VARS == [2, 3, 3]
    path = write_vcf(WORKED_EXAMPLE_LINES)

    at_three = calc_pairwise_kosman_dists(open_vcf(path), min_num_snps=3)
    at_four = calc_pairwise_kosman_dists(open_vcf(path), min_num_snps=4)

    assert math.isnan(at_three.dist_vector[0])
    assert list(at_three.dist_vector[1:]) == WORKED_EXAMPLE_DISTS[1:]
    assert numpy.all(numpy.isnan(at_four.dist_vector))


def test_the_distances_of_two_individuals_are_the_ones_pynei_asserts(
    vcf_of,
) -> None:
    """The genotypes of pyNei's `test_kosman_2_indis`, its four individuals
    in one dataset.

    pyNei makes the three checks one pair at a time, at a class that popnei
    has nothing like, and popnei makes them at the public function: the
    pairs of a dataset are worked out from their own genotypes, so the three
    distances are the ones of its three datasets of two individuals.
    """
    gts = _of_individuals(
        KOSMAN_2_INDIS_A, KOSMAN_2_INDIS_B, KOSMAN_2_INDIS_C, KOSMAN_2_INDIS_D
    )

    dists = calc_pairwise_kosman_dists(
        open_vcf(vcf_of(gts, ["a", "b", "c", "d"]))
    ).square_dists

    assert dists.loc["a", "b"] == DIST_A_B
    assert dists.loc["c", "d"] == DIST_C_D
    assert dists.loc["b", "d"] == DIST_B_D


def test_a_variant_missing_in_both_individuals_counts_as_one_missing_in_one(
    vcf_of,
) -> None:
    """The two datasets of pyNei's `test_kosman_missing`.

    A variant counts for a pair only when both genotypes are called, so
    making the genotype of the second individual missing at the variant
    where the first one is already missing changes nothing.
    """
    in_each = _of_individuals(*MISSING_IN_EACH)
    in_both = _of_individuals(*MISSING_IN_BOTH)

    of_each = calc_pairwise_kosman_dists(open_vcf(vcf_of(in_each, name="each.vcf")))
    of_both = calc_pairwise_kosman_dists(open_vcf(vcf_of(in_both, name="both.vcf")))

    assert of_each.dist_vector[0] == DIST_A_B
    assert of_each.dist_vector[0] == of_both.dist_vector[0]


def test_the_pairwise_distances_are_the_vector_pynei_asserts(vcf_of) -> None:
    """The four individuals of pyNei's `test_kosman_pairwise`, one of whose
    genotypes holds the alleles 1 and 2.

    The vector is in the order of the upper triangle of the square matrix,
    (a, b), (a, c), (a, d), (b, c), (b, d), (c, d).
    """
    gts = _of_individuals(
        KOSMAN_2_INDIS_A, KOSMAN_PAIRWISE_B, KOSMAN_2_INDIS_C, KOSMAN_2_INDIS_D
    )

    dists = calc_pairwise_kosman_dists(open_vcf(vcf_of(gts, ["a", "b", "c", "d"])))

    assert list(dists.dist_vector) == KOSMAN_PAIRWISE_DISTS


def test_a_filter_that_takes_nothing_out_leaves_the_distances_and_counts_all(
    vcf_of,
) -> None:
    """The dataset of pyNei's `test_kosman_pairwise_with_filtered_vars`, 60
    variants of 30 individuals with every genotype called, with a missing
    data filter at 1 as a step of the `Variants`.

    The filter takes nothing out, so the distances are those of the
    `Variants` without it, twice in a row: the second call makes a pass of
    its own over the source. The counts of each pass say that the filter was
    given every variant and kept every one of them.
    """
    rng = numpy.random.default_rng(FILTERED_VARS_SEED)
    gts = rng.integers(
        0, 2, size=(FILTERED_VARS_NUM_VARS, FILTERED_VARS_NUM_INDIVIDUALS, 2)
    )
    path = vcf_of(gts)
    expected = calc_pairwise_kosman_dists(open_vcf(path)).dist_vector

    filtered = open_vcf(path)
    filtered.filter_by_missing_data(1)

    for _ in range(2):
        dists = calc_pairwise_kosman_dists(filtered)
        assert numpy.array_equal(dists.dist_vector, expected)
        assert dists.pass_stats == PassStats(
            num_vars=FILTERED_VARS_NUM_VARS,
            filtering={
                "missing_data": FilteringStats(
                    vars_processed=FILTERED_VARS_NUM_VARS,
                    vars_kept=FILTERED_VARS_NUM_VARS,
                )
            },
        )


def test_the_distances_are_over_the_variants_the_steps_kept(write_vcf) -> None:
    """Five variants of three individuals, two of which have a missing
    genotype, and a filter that keeps the variants with none.

    The distances are those of the three variants that passed, and the
    counts say so: `num_vars` is the 3 the calculation took and not the 5
    the source gave, which the filter reports beside it.
    """
    variants = open_vcf(
        write_vcf(
            [
                "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
                "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t0/1\t1/1",
                "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1",
                "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t0/0",
                "chr1\t50\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0",
            ]
        )
    )
    variants.filter_by_missing_data(0)

    dists = calc_pairwise_kosman_dists(variants)

    # Over the three variants that passed, the two individuals of each pair
    # share one allele at the first and the last and hold the same genotype
    # at the middle one, but for ind1 and ind3, which share no allele at the
    # first and the last: 2 over 2 x 3, 4 over 2 x 3 and 2 over 2 x 3.
    assert list(dists.dist_vector) == [1 / 3, 2 / 3, 1 / 3]
    assert dists.pass_stats == PassStats(
        num_vars=3,
        filtering={"missing_data": FilteringStats(vars_processed=5, vars_kept=3)},
    )


def test_a_pair_called_together_at_no_variant_has_no_distance(vcf_of) -> None:
    """Three individuals over two variants, where the first two are never
    called at the same variant.

    Nothing else is affected by such a pair: the other pairs of the two
    individuals keep their distances, the square matrix has NaN in the two
    cells of the pair, and its diagonal is 0.
    """
    gts = [
        [[-1, -1], [0, 0], [0, 1]],
        [[0, 0], [-1, -1], [1, 1]],
    ]

    dists = calc_pairwise_kosman_dists(open_vcf(vcf_of(gts, ["a", "b", "c"])))

    assert math.isnan(dists.dist_vector[0])
    # a and c share no allele at the variant both are called at, and b and c
    # share one of the two.
    assert list(dists.dist_vector[1:]) == [1.0, 0.5]
    square = dists.square_dists
    assert math.isnan(square.loc["a", "b"])
    assert math.isnan(square.loc["b", "a"])
    assert list(numpy.diag(square.values)) == [0.0, 0.0, 0.0]


def test_one_individual_gives_an_empty_vector(vcf_of) -> None:
    """A dataset of one individual, which makes no pair.

    Its square matrix is the one cell of its diagonal, which is 0, and what
    a user prints of the result counts that individual in the singular.
    """
    gts = [[[0, 0]], [[0, 1]]]

    dists = calc_pairwise_kosman_dists(open_vcf(vcf_of(gts, ["only"])))

    assert dists.dist_vector.shape == (0,)
    assert dists.names == ("only",)
    assert dists.square_dists.values.tolist() == [[0.0]]
    assert dists.triang_list_of_lists == [[0.0]]
    assert repr(dists) == (
        "<Distances of 1 individual, 0 pairs, with the counts of its pass>"
    )
    assert repr(Distances([], names=["only"])) == (
        "<Distances of 1 individual, 0 pairs>"
    )


def test_the_counts_of_the_pass_hold_the_variants_of_the_dataset() -> None:
    """The panel, whose 1200 variants the calculation was given.

    No filter is on the `Variants`, so there is no count of a filter to
    read, and `num_vars` is what the calculation took.
    """
    name, ploidy, num_vars, _ = PANEL

    dists = calc_pairwise_kosman_dists(open_vcf(_reference(name), ploidy=ploidy))

    assert dists.pass_stats == PassStats(num_vars=num_vars, filtering={})


def test_the_variants_are_as_they_were_after_the_calculation(
    write_vcf,
) -> None:
    """One pass over the source, and a `Variants` that nothing was added to.

    A user gives one `Variants` to one calculation after another, so the
    steps have to be what they were and the second call has to give what the
    first gave.
    """
    variants = open_vcf(write_vcf(WORKED_EXAMPLE_LINES))
    variants.filter_by_missing_data(1)
    steps_before = variants.steps

    first = calc_pairwise_kosman_dists(variants)
    second = calc_pairwise_kosman_dists(variants)

    assert variants.steps == steps_before
    assert len(variants.steps) == 1
    assert numpy.array_equal(first.dist_vector, second.dist_vector)


def test_a_negative_min_num_snps_is_refused(write_vcf) -> None:
    """A `min_num_snps` below 0, which pyNei takes and does nothing with.

    How many variants a pair needs is a count, so a negative one says
    nothing a user can mean, and popnei says so at the call. So do 2.5
    variants and `True`, which Python would pass on as the number 1, and a
    number of variants above the 4294967295 popnei counts
    for a pair, which no pair could reach and which would otherwise give a
    vector of nothing but NaN.
    """
    variants = open_vcf(write_vcf(WORKED_EXAMPLE_LINES))

    with pytest.raises(ValueError, match="min_num_snps"):
        calc_pairwise_kosman_dists(variants, min_num_snps=-1)
    with pytest.raises(TypeError, match="min_num_snps"):
        calc_pairwise_kosman_dists(variants, min_num_snps=2.5)
    with pytest.raises(TypeError, match="min_num_snps"):
        calc_pairwise_kosman_dists(variants, min_num_snps=True)
    with pytest.raises(ValueError, match="4294967295"):
        calc_pairwise_kosman_dists(variants, min_num_snps=2**40)


def _what_it_said_of(refusal: pytest.ExceptionInfo, path: Path) -> str:
    """What a refusal says after the path of the file, which every message
    of a file starts with."""
    message = str(refusal.value)
    assert message.startswith(f"{path}: ")
    return message.removeprefix(f"{path}: ")


def test_a_source_with_no_variant_is_refused(write_vcf) -> None:
    """A VCF whose header names three individuals and that has no data line.

    A calculation over no variant gives no number, so it is a wrong input
    and not a result of NaN, and the message says which of the two happened,
    the source having none or the steps keeping none. This one is the first,
    and the whole sentence is asserted, because every word of it is what
    tells the two apart.
    """
    path = write_vcf([])

    with pytest.raises(ValueError) as refusal:
        calc_pairwise_kosman_dists(open_vcf(path))

    assert _what_it_said_of(refusal, path) == (
        "the source has no variant, and a calculation needs 1 variant at least"
    )


def test_a_source_with_no_variant_is_told_apart_from_steps_that_kept_none(
    write_vcf,
) -> None:
    """The same VCF with no data line, with a filter on it.

    The filter was given no variant, which is what says that the source is
    the one that had none: the counts come after the sentence that says it.
    """
    path = write_vcf([])
    variants = open_vcf(path)
    variants.filter_by_missing_data(0)

    with pytest.raises(ValueError) as refusal:
        calc_pairwise_kosman_dists(variants)

    assert _what_it_said_of(refusal, path) == (
        "the source has no variant, and a calculation needs 1 variant at "
        "least: the filter `missing_data` was given 0 variants and kept 0"
    )


def test_steps_that_kept_no_variant_are_refused_with_what_each_filter_counted(
    write_vcf,
) -> None:
    """Four variants, each with one genotype of the three missing, and a
    missing data filter that keeps the variants with no missing genotype,
    with a second filter after it.

    The counts of a pass that could not finish are otherwise lost, so the
    message carries them, in the order of the steps: the first filter was
    given the four variants of the source and kept none, and the second was
    given none. The number of variants the source gave is in the sentence
    itself, which is what says that the steps and not the source are what
    left the calculation with nothing.
    """
    path = write_vcf(
        [
            "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t./.",
            "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
            "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t./.\t0/1\t1/1",
            "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t1/1",
        ]
    )
    variants = open_vcf(path)
    variants.filter_by_missing_data(0)

    with pytest.raises(ValueError) as refusal:
        calc_pairwise_kosman_dists(variants)
    assert _what_it_said_of(refusal, path) == (
        "the steps kept no variant of the 4 the source gave, and a "
        "calculation needs 1 variant at least: the filter `missing_data` was "
        "given 4 variants and kept 0"
    )

    with_two = open_vcf(path)
    with_two.filter_by_missing_data(0)
    with_two.filter_by_maf(0.5)

    with pytest.raises(ValueError) as refusal:
        calc_pairwise_kosman_dists(with_two)
    assert _what_it_said_of(refusal, path) == (
        "the steps kept no variant of the 4 the source gave, and a "
        "calculation needs 1 variant at least: the filter `missing_data` was "
        "given 4 variants and kept 0, the filter `maf` was given 0 variants "
        "and kept 0"
    )


def test_a_vector_of_no_number_of_individuals_is_refused() -> None:
    """Four distances, which are the pairs of no number of individuals.

    pyNei takes them and builds the three individuals of the three first
    values, with the fourth in no cell of its square matrix. A square matrix
    given in the place of the vector is refused too, and names the way to
    build a result from one.
    """
    with pytest.raises(ValueError) as refusal:
        Distances(VECTOR_OF_NO_NUMBER_OF_INDIVIDUALS)
    # The numbers of pairs around 4 are those of 3 and of 4 individuals.
    assert "3 for 3 individuals" in str(refusal.value)
    assert "6 for 4" in str(refusal.value)
    with pytest.raises(ValueError, match="from_square_dists"):
        Distances([[0.0, 0.5], [0.5, 0.0]])


def test_what_a_distances_cannot_be_built_from_is_refused_by_its_argument() -> None:
    """Distances that are no numbers, names that are no sequence of names,
    and a square matrix that is no frame.

    Each of the three used to come out as the message of numpy or of Python,
    "could not convert string to float", "'int' object is not iterable" and
    an `AttributeError`, none of which names the argument the user wrote.
    """
    with pytest.raises(ValueError, match="dist_vector"):
        Distances(numpy.array(["a", "b", "c"]))
    with pytest.raises(TypeError, match="names"):
        Distances([0.5, 0.75, 0.3], names=3)
    with pytest.raises(TypeError, match="square_dists"):
        Distances.from_square_dists(numpy.zeros((3, 3)))


def test_the_triangle_of_a_distances_is_the_lower_one_of_its_square_matrix() -> None:
    """The four individuals of "What pyNei does that is odd", each distance
    saying which pair it is.

    pyNei cuts the vector into runs of 1, 2, 3 values, which are not the
    rows of the lower triangle from four individuals on: it gives the row of
    the third individual as 2, 3 where the square matrix has 2, 12.
    """
    dists = Distances(ODD_VECTOR, names=["a", "b", "c", "d"])

    assert dists.triang_list_of_lists == ODD_TRIANGLE
    # Every value of it is a distance, the 0 of the diagonal included, so a
    # user who writes it to a file of their own gets one kind of number.
    assert all(
        isinstance(dist, float) for row in dists.triang_list_of_lists for dist in row
    )
    square = dists.square_dists
    assert [square.loc["c", "a"], square.loc["c", "b"]] == [2, 12]


def test_a_distances_is_built_from_a_square_matrix_and_gives_it_back() -> None:
    """A square frame of distances calculated elsewhere, and the vector of
    the upper triangle popnei keeps them in.

    The names come from the index of the frame, and a frame whose two sides
    are of different sizes is refused.
    """
    names = ["a", "b", "c", "d"]
    dists = Distances(ODD_VECTOR, names=names)

    again = Distances.from_square_dists(dists.square_dists)

    assert list(again.dist_vector) == ODD_VECTOR
    assert again.names == tuple(names)
    assert again.pass_stats is None
    assert again.square_dists.equals(dists.square_dists)
    with pytest.raises(ValueError, match="square"):
        Distances.from_square_dists(pandas.DataFrame([[0.0, 0.1, 0.2]]))


def test_a_square_matrix_whose_two_sides_are_not_in_one_order_is_refused() -> None:
    """A frame whose rows were sorted and whose columns were not, which
    `square.loc[["c", "b", "a"]]` gives.

    Its cells are no longer the distance of the pair of their row and their
    column, and reading the upper triangle of it gave a wrong vector in
    silence: for the three individuals a, b, c at 0.1, 0.2 and 0.3 it gave
    0.3, 0.0, 0.3 under the names c, b, a, where the distance of c and a is
    0.2 and no pair has a distance of 0.
    """
    dists = Distances([0.1, 0.2, 0.3], names=["a", "b", "c"])
    turned_around = dists.square_dists.loc[["c", "b", "a"]]

    with pytest.raises(ValueError) as refusal:
        Distances.from_square_dists(turned_around)

    message = str(refusal.value)
    assert "'c'" in message
    assert "'a'" in message
    assert "0" in message
    # The frame whose two sides were sorted together is taken, and its
    # distances are those of the pairs of its own order: c with b, c with a
    # and b with a.
    sorted_together = dists.square_dists.loc[["c", "b", "a"], ["c", "b", "a"]]
    of_the_other_order = Distances.from_square_dists(sorted_together)
    assert of_the_other_order.names == ("c", "b", "a")
    assert list(of_the_other_order.dist_vector) == [0.3, 0.2, 0.1]


def test_a_distances_built_with_no_name_names_its_individuals_by_their_place() -> None:
    """Three distances and no name, which pyNei names 0, 1 and 2 too.

    The names are a tuple, as the names of individuals are everywhere in
    popnei, where pyNei has a numpy array of them.
    """
    dists = Distances([0.5, 0.75, 0.3])

    assert dists.names == (0, 1, 2)
    assert isinstance(dists.names, tuple)
    assert list(dists.square_dists.index) == [0, 1, 2]
    with pytest.raises(ValueError) as refusal:
        Distances([0.5, 0.75, 0.3], names=["a", "b"])
    assert "2 names" in str(refusal.value)
    assert "3 individuals" in str(refusal.value)


def test_an_empty_vector_is_of_no_individual_or_of_one() -> None:
    """No distance, which is what one individual gives and what no
    individual would give.

    The two cannot be told apart from the vector, so both are taken, and a
    number of names that is neither says so.
    """
    assert Distances([], names=[]).names == ()
    assert Distances([], names=["only"]).names == ("only",)
    assert Distances([]).names == (0,)

    with pytest.raises(ValueError) as refusal:
        Distances([], names=["a", "b"])
    assert "0 or of 1 individual" in str(refusal.value)


def test_the_vector_of_a_distances_cannot_be_written_into(
    write_vcf,
) -> None:
    """The array of a result, which nothing changes after it was given.

    A user who keeps a result and works from it reads the same numbers as
    long as they hold it.
    """
    dists = calc_pairwise_kosman_dists(open_vcf(write_vcf(WORKED_EXAMPLE_LINES)))

    assert dists.dist_vector.dtype == numpy.float64
    with pytest.raises(ValueError):
        dists.dist_vector[0] = 0.0


def test_a_haploid_and_a_tetraploid_variants_are_taken(vcf_of) -> None:
    """A tetraploid and a haploid `Variants`, which pyNei refuses with "Only
    diploid are allowed".

    Every ploidy is taken, and the ploidy of the source is what the sum of d
    is divided by: two haploid individuals that hold a different allele at
    every variant are 1 apart, and two tetraploid ones that share two of
    their four alleles are 0.5 apart.
    """
    haploid = [[[0], [1]], [[0], [1]]]
    tetraploid = [[[0, 0, 0, 1], [0, 1, 1, 1]], [[0, 0, 0, 1], [0, 1, 1, 1]]]

    of_haploid = calc_pairwise_kosman_dists(
        open_vcf(vcf_of(haploid, ["h0", "h1"], name="h.vcf"), ploidy=1)
    )
    of_tetraploid = calc_pairwise_kosman_dists(
        open_vcf(vcf_of(tetraploid, ["t0", "t1"], name="t.vcf"), ploidy=4)
    )

    assert list(of_haploid.dist_vector) == [1.0]
    assert list(of_tetraploid.dist_vector) == [0.5]


def test_a_distances_holds_the_names_as_a_tuple_and_the_counts_of_its_pass(
    write_vcf,
) -> None:
    """What a result carries besides the distances.

    The names are a tuple, so what a user reads is not a list of popnei's
    that they could write into, and `pass_stats` is `None` in a `Distances`
    that was built from distances calculated elsewhere.
    """
    dists = calc_pairwise_kosman_dists(open_vcf(write_vcf(WORKED_EXAMPLE_LINES)))

    assert isinstance(dists.names, tuple)
    assert dists.names == ("ind1", "ind2", "ind3")
    assert dists.pass_stats == PassStats(num_vars=4, filtering={})
    assert Distances(ODD_VECTOR).pass_stats is None


def test_the_square_matrix_is_indexed_by_the_names_on_both_sides(
    write_vcf,
) -> None:
    """The N x N frame a user gives to a tree or a principal coordinate
    analysis.

    It is symmetrical, its diagonal is 0, and both of its indexes are the
    names of the individuals in the order the source has them.
    """
    dists = calc_pairwise_kosman_dists(open_vcf(write_vcf(WORKED_EXAMPLE_LINES)))

    square = dists.square_dists

    assert list(square.index) == ["ind1", "ind2", "ind3"]
    assert list(square.columns) == ["ind1", "ind2", "ind3"]
    assert square.loc["ind1", "ind2"] == WORKED_EXAMPLE_DISTS[0]
    assert square.loc["ind2", "ind1"] == WORKED_EXAMPLE_DISTS[0]
    assert list(numpy.diag(square.values)) == [0.0, 0.0, 0.0]


def test_a_variants_of_another_source_gives_the_same_distances(
    tmp_path: Path, write_vcf
) -> None:
    """The worked example read from a VCF and from the vars file written
    from it.

    The calculation reads whatever source the `Variants` holds, and the two
    integers of every pair are the same whichever file the genotypes came
    from and however the blocks fall.
    """
    path = tmp_path / "worked.vars"
    write_vars(open_vcf(write_vcf(WORKED_EXAMPLE_LINES)), path, 2)

    of_the_vars_file = calc_pairwise_kosman_dists(open_vars(path))

    assert list(of_the_vars_file.dist_vector) == WORKED_EXAMPLE_DISTS
    assert of_the_vars_file.names == ("ind1", "ind2", "ind3")


def test_a_variants_is_what_the_calculation_takes(write_vcf) -> None:
    """A `Variants` and not a path: the mistake that is easiest to make."""
    with pytest.raises(TypeError, match="open_vcf"):
        calc_pairwise_kosman_dists(str(write_vcf(WORKED_EXAMPLE_LINES)))


def test_the_result_is_of_the_individuals_of_the_source(write_vcf) -> None:
    """The names of a result and those of the `Variants` it came from."""
    variants: Variants = open_vcf(write_vcf(WORKED_EXAMPLE_LINES))

    dists = calc_pairwise_kosman_dists(variants)

    assert dists.names == variants.individuals
