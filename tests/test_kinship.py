"""The kinship of every pair of individuals, from Python.

`docs/specs/kinship.md` has the matrix, the `Kinship` it comes in and the
numbers this file asserts. Every literal here is the spec's or is read from
`tests/reference/kinship/`, which
`tests/reference/kinship/make_reference.py` writes: the panel with every
genotype called as a gzipped VCF, and, for it and for the panel of
`tests/reference/dists/` with 3 in 100 of its genotypes missing, the matrix
that `plink2 --make-rel square` wrote for it with the individuals of that
matrix in its order.

The two checks the spec asks for are against plink2 v2.0.0-a.7.7, within
1e-5 absolute, which is one unit of the last of the six significant digits
plink2 prints for an entry near 1; and against pyNei at commit ef0ca6e,
which `pyproject.toml` names, within 1e-12 relative, because the two
libraries add the variants of a pair in different orders. pyNei reads the
same VCFs, so both libraries are given the same genotypes.

The datasets that no file holds are written here as small VCFs: the worked
example of the spec, the pair of individuals that were never called
together, and the variant of three alleles.
"""

import gzip
from collections.abc import Sequence
from pathlib import Path

import numpy
import pandas
import pytest
from popnei import Kinship, calc_kinship, open_vcf
from pynei import vars_from_vcf
from pynei.gwas import calc_kinship as pynei_kinship

REFERENCE_KINSHIP_DIR = Path(__file__).parent / "reference" / "kinship"
REFERENCE_DISTS_DIR = Path(__file__).parent / "reference" / "dists"

# The two panels of "How it is verified" of the spec, the same 200
# individuals and 1200 biallelic diploid variants twice: `panel_called` with
# every genotype called, and `panel` with 3 in 100 of them missing whole,
# which is the one that makes the denominator of a pair differ from the
# denominator of another.
PANELS = ("panel_called", "panel")
PANEL_NUM_VARS = 1200
PANEL_NUM_INDIVIDUALS = 200

# What the tests compare within, and why. The whole of each matrix is held
# to 1e-12 relative against the float64 plink2 wrote, which "How it is
# verified" of the spec gives: the worst entry of the two panels is 3.3e-13
# away as a ratio, at the smallest entries, where the two libraries add the
# same products in a different order, so 1e-12 is two to three times the
# worst and not a wide margin. pyNei is held to the same bound for the same
# reason. The one entry read from the text plink2 prints is held to 1e-5
# absolute, which is what six significant digits of an entry near 1 allow.
OF_PLINK2 = 1e-12
OF_PYNEI = 1e-12
OF_THE_PRINTED_DIGITS = 1e-5

# The individuals 10 to 49 of the panel with every genotype called, the 40
# of "Its Python function, and its TypeScript one", and what their kinship
# gives: 1195 of the 1200 variants have variance among them, and their
# matrix is as far as 0.129 from the same 40 rows and columns of the kinship
# of all 200.
THE_FORTY = slice(10, 50)
NUM_VARS_OF_THE_FORTY = 1195
FURTHEST_FROM_THE_SLICE = 0.129

# The entry of the pair `s000`, `s001`, two full sibs of the panel with every
# genotype called, from the table of the spec, which plink2 wrote. It is
# looked up by name, so it says that the names index the rows they belong to.
OF_TWO_FULL_SIBS = 0.648081

# The worked example of "How it is verified": 4 diploid individuals and 4
# variants, of which the third has no variance, every individual being
# heterozygous, and the fourth has one allele. The genotype of `i2` at the
# second variant is missing, which is what leaves the pairs of `i2` with one
# variant where the others have two.
WORKED_EXAMPLE_GTS = [
    [[0, 0], [0, 1], [1, 1], [0, 1]],
    [[0, 0], [0, 1], [-1, -1], [1, 1]],
    [[0, 1], [0, 1], [0, 1], [0, 1]],
    [[0, 0], [0, 0], [0, 0], [0, 0]],
]
WORKED_EXAMPLE_INDIVIDUALS = ("i0", "i1", "i2", "i3")
WORKED_EXAMPLE_NUM_VARS = 2
WORKED_EXAMPLE_MATRIX = [
    [2.0, 0.0, -2.0, -1.0],
    [0.0, 0.0, 0.0, 0.0],
    [-2.0, 0.0, 2.0, 0.0],
    [-1.0, 0.0, 0.0, 1.0],
]
# The same four variants over `i0` and `i3` alone: both variants have
# variance among those two, every denominator is 2, and the matrix is 4/3 on
# the diagonal and -4/3 off it, which is not the two rows and columns of the
# matrix above.
WORKED_EXAMPLE_OF_TWO = [[4 / 3, -4 / 3], [-4 / 3, 4 / 3]]
# pyNei gives every entry of both as a whole number within 4.4e-16, so the
# spec asserts them within 1e-12 absolute.
OF_THE_WORKED_EXAMPLE = 1e-12

# Two variants of three individuals in which the first and the third are
# never called together, so the sum of that pair would be divided by no
# variant at all. Each variant has variance among the individuals whose
# genotype is called, so neither is left out before the pair is reached.
NEVER_CALLED_TOGETHER_GTS = [
    [[0, 0], [1, 1], [-1, -1]],
    [[-1, -1], [0, 0], [1, 1]],
]
# Their names, which tell a message that names the individuals from one that
# names the places they are at: `first` is at the place 0 and `third` at 2.
NEVER_CALLED_TOGETHER_INDIVIDUALS = ("first", "second", "third")

# Two variants of three individuals of which the second has three alleles
# among its called genotypes, `0/1`, `1/2` and `2/2`.
THREE_ALLELES_GTS = [
    [[0, 0], [0, 1], [1, 1]],
    [[0, 1], [1, 2], [2, 2]],
]

# One more individual than the 46340 a kinship is taken of, whose square is
# the 2147483647 values the linear algebra counts in.
TOO_MANY_INDIVIDUALS = 46341

_NUCLEOTIDES = "ACGT"


def _vcf_text(gts: numpy.ndarray, individuals: Sequence[str]) -> str:
    """The genotypes `gts`, variants x individuals x ploidy, as a VCF.

    Every variant declares the alleles of the whole dataset, one nucleotide
    each, and an allele below 0 is written as a dot: which letter an allele
    gets changes no entry of the kinship, which is built from the dosages,
    and REF and ALT only name them.
    """
    alleles = _NUCLEOTIDES[: max(2, int(gts.max()) + 1)]
    lines = [
        "##fileformat=VCFv4.4",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + list(individuals)
        ),
    ]
    for index in range(gts.shape[0]):
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


def _panel(name: str) -> Path:
    """The gzipped VCF of one of the two panels. The one with genotypes
    missing is the panel of `docs/specs/dists.md`, which was in the
    repository before the kinship and is not written twice."""
    if name == "panel":
        return REFERENCE_DISTS_DIR / "panel.vcf.gz"
    return REFERENCE_KINSHIP_DIR / f"{name}.vcf.gz"


def _plink2_matrix(name: str) -> numpy.ndarray:
    """The 40000 entries of a panel as plink2 holds them, the little endian
    float64 of `--make-rel square bin`, row after row.

    The `<name>.plink2.rel.gz` beside it is the same matrix as text of six
    significant digits, which is what the table of literals of the spec is
    read from and what one entry of it is asserted against here. The whole
    of the matrix is compared with these bits: six digits round an entry
    near 1 by up to 5e-6, so a comparison with the text within 1e-5 passes
    an error of up to 5e-6 and has no room left to find one in.
    """
    with gzip.open(REFERENCE_KINSHIP_DIR / f"{name}.plink2.rel.bin.gz", "rb") as stored:
        entries = numpy.frombuffer(stored.read(), dtype="<f8")
    return entries.reshape(PANEL_NUM_INDIVIDUALS, PANEL_NUM_INDIVIDUALS)


def _plink2_individuals(name: str) -> tuple[str, ...]:
    """The individuals of that matrix in its order, which plink2 wrote beside
    it under `#IID`."""
    lines = (REFERENCE_KINSHIP_DIR / f"{name}.plink2.rel.id").read_text().splitlines()
    return tuple(line for line in lines[1:] if line)


@pytest.mark.parametrize("name", PANELS)
def test_every_entry_of_a_panel_is_the_one_plink2_wrote(name: str) -> None:
    """The 40000 entries of each panel against the bits of plink2's matrix.

    The individuals are compared with the ones plink2 wrote beside its
    matrix, so that the entries are matched pair by pair and not by their
    place alone, and the entry of the two full sibs `s000` and `s001` is
    looked up by name in the digits the spec's table prints.
    """
    kinship = calc_kinship(open_vcf(_panel(name)))
    of_plink2 = _plink2_matrix(name)

    assert kinship.num_vars == PANEL_NUM_VARS
    assert kinship.individuals == _plink2_individuals(name)
    assert kinship.matrix.shape == (PANEL_NUM_INDIVIDUALS, PANEL_NUM_INDIVIDUALS)
    assert of_plink2.shape == kinship.matrix.shape
    numpy.testing.assert_allclose(
        kinship.matrix.to_numpy(), of_plink2, rtol=OF_PLINK2, atol=0
    )
    if name == "panel_called":
        assert kinship.matrix.loc["s000", "s001"] == pytest.approx(
            OF_TWO_FULL_SIBS, abs=OF_THE_PRINTED_DIGITS
        )


@pytest.mark.parametrize("name", PANELS)
def test_every_entry_of_a_panel_is_pyneis(name: str) -> None:
    """Both libraries on the same VCF, entry by entry and variant count.

    pyNei gives no counts of the pass, so what is compared is the matrix, the
    individuals and how many variants each library used.
    """
    ours = calc_kinship(open_vcf(_panel(name)))
    theirs = pynei_kinship(vars_from_vcf(_panel(name)))

    assert ours.num_vars == theirs.num_vars
    assert ours.individuals == tuple(theirs.samples)
    numpy.testing.assert_allclose(
        ours.matrix.to_numpy(), theirs.matrix.to_numpy(), rtol=OF_PYNEI, atol=0
    )


def test_the_kinship_of_40_individuals_is_not_the_40_rows_of_the_whole_panel() -> None:
    """`individuals` takes the frequencies of those individuals.

    The 40 individuals 10 to 49 of the panel with every genotype called have
    1195 of its 1200 variants with variance among them, and their matrix is
    as far as 0.129 from the same 40 rows and columns of the kinship of all
    200, which is what the frequencies, the means and the denominators being
    theirs and not the panel's does.
    """
    whole = calc_kinship(open_vcf(_panel("panel_called")))
    the_forty = whole.individuals[THE_FORTY]
    of_the_forty = calc_kinship(open_vcf(_panel("panel_called")), individuals=the_forty)

    assert whole.num_vars == PANEL_NUM_VARS
    assert of_the_forty.num_vars == NUM_VARS_OF_THE_FORTY
    assert of_the_forty.individuals == the_forty
    sliced = whole.matrix.loc[list(the_forty), list(the_forty)].to_numpy()
    furthest = numpy.abs(of_the_forty.matrix.to_numpy() - sliced).max()
    assert furthest == pytest.approx(FURTHEST_FROM_THE_SLICE, abs=5e-4)


def test_the_kinship_of_40_individuals_is_pyneis() -> None:
    """The same 40 individuals in both libraries, which pyNei calculates the
    frequencies of over those 40 as well."""
    ours = calc_kinship(open_vcf(_panel("panel_called")))
    the_forty = list(ours.individuals[THE_FORTY])
    ours = calc_kinship(open_vcf(_panel("panel_called")), individuals=the_forty)
    theirs = pynei_kinship(vars_from_vcf(_panel("panel_called")), samples=the_forty)

    assert ours.num_vars == theirs.num_vars
    assert ours.individuals == tuple(theirs.samples)
    numpy.testing.assert_allclose(
        ours.matrix.to_numpy(), theirs.matrix.to_numpy(), rtol=OF_PYNEI, atol=0
    )


def test_the_individuals_are_taken_in_the_order_they_were_named(vcf_of) -> None:
    """`individuals` in another order than the file's gives the matrix in
    that order, and every entry is the entry of its pair."""
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))
    turned = calc_kinship(variants, individuals=("i3", "i2", "i1", "i0"))
    as_they_come = calc_kinship(
        open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))
    )

    assert turned.individuals == ("i3", "i2", "i1", "i0")
    assert turned.num_vars == as_they_come.num_vars
    for one in WORKED_EXAMPLE_INDIVIDUALS:
        for other in WORKED_EXAMPLE_INDIVIDUALS:
            assert turned.matrix.loc[one, other] == pytest.approx(
                as_they_come.matrix.loc[one, other], abs=OF_THE_WORKED_EXAMPLE
            )


def test_the_worked_example_gives_the_matrix_of_whole_numbers(vcf_of) -> None:
    """The four variants of four individuals of "How it is verified".

    Two of the four variants have variance, so `num_vars` is 2 where the pass
    gave 4; the pairs of `i2`, whose genotype is missing at the second
    variant, are divided by 1 and the others by 2; and `i1` is 0 against
    everyone, both of its genotypes being heterozygous and standardized to 0.
    """
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))
    kinship = calc_kinship(variants)

    assert kinship.individuals == WORKED_EXAMPLE_INDIVIDUALS
    assert kinship.num_vars == WORKED_EXAMPLE_NUM_VARS
    assert kinship.pass_stats.num_vars == len(WORKED_EXAMPLE_GTS)
    assert kinship.pass_stats.filtering == {}
    numpy.testing.assert_allclose(
        kinship.matrix.to_numpy(),
        numpy.array(WORKED_EXAMPLE_MATRIX),
        rtol=0,
        atol=OF_THE_WORKED_EXAMPLE,
    )


def test_the_worked_example_of_two_individuals_has_their_frequencies(vcf_of) -> None:
    """The same four variants over `i0` and `i3` alone.

    The first variant has the dosages 0 and 1 among those two, so its mean
    dosage is 0.5 where over the four individuals it is 1, and both variants
    are kept: the matrix is 4/3 on the diagonal and -4/3 off it, and it is
    not the two rows and columns of the matrix of the four.
    """
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))
    kinship = calc_kinship(variants, individuals=("i0", "i3"))

    assert kinship.individuals == ("i0", "i3")
    assert kinship.num_vars == WORKED_EXAMPLE_NUM_VARS
    numpy.testing.assert_allclose(
        kinship.matrix.to_numpy(),
        numpy.array(WORKED_EXAMPLE_OF_TWO),
        rtol=0,
        atol=OF_THE_WORKED_EXAMPLE,
    )


def test_the_counts_of_the_pass_are_those_of_the_filter_and_not_of_the_kinship(
    vcf_of,
) -> None:
    """The worked example with the missing data filter at 0 on it.

    The filter takes the second variant away, so the pass gives 3 variants
    and the kinship is left with the first alone, the other two having no
    variance: the counts of the pass hold what the filter was given and
    kept, 4 and 3, and `num_vars` is the 1 variant that was used.
    """
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))
    variants.filter_by_missing_data(0)
    kinship = calc_kinship(variants)

    assert kinship.num_vars == 1
    assert kinship.pass_stats.num_vars == 3
    assert kinship.pass_stats.filtering["missing_data"].vars_processed == 4
    assert kinship.pass_stats.filtering["missing_data"].vars_kept == 3
    # The first variant alone, whose standardized dosages are -sqrt(2), 0,
    # sqrt(2) and 0, each pair divided by the one variant it was called at.
    numpy.testing.assert_allclose(
        kinship.matrix.to_numpy(),
        numpy.array(
            [
                [2.0, 0.0, -2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
                [-2.0, 0.0, 2.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ]
        ),
        rtol=0,
        atol=OF_THE_WORKED_EXAMPLE,
    )


def test_a_variant_of_three_alleles_is_refused_and_transformed_when_asked(
    vcf_of,
) -> None:
    """A variant with three alleles among its called genotypes.

    Without `transform_to_biallelic` it is a `ValueError` that says which
    variant it is and which argument to pass, and with it every allele that
    is not the major one counts the same, which is what pyNei does to every
    variant silently: the two libraries then give the same matrix.
    """
    path = vcf_of(THREE_ALLELES_GTS)
    with pytest.raises(ValueError, match="transform_to_biallelic"):
        calc_kinship(open_vcf(path))

    ours = calc_kinship(open_vcf(path), transform_to_biallelic=True)
    theirs = pynei_kinship(vars_from_vcf(path))

    assert ours.num_vars == theirs.num_vars
    numpy.testing.assert_allclose(
        ours.matrix.to_numpy(), theirs.matrix.to_numpy(), rtol=OF_PYNEI, atol=0
    )


def test_a_pair_with_no_variant_called_in_both_names_the_two(vcf_of) -> None:
    """Two variants of three individuals in which the first and the third are
    never called together.

    Their entry would be divided by no variant at all. popnei raises a
    `ValueError` naming the two individuals, and not the places they are at,
    and how many variants each of them has called, where pyNei divides and
    leaves a NaN in the matrix. The file is named as it is in every error of
    one.
    """
    path = vcf_of(NEVER_CALLED_TOGETHER_GTS, NEVER_CALLED_TOGETHER_INDIVIDUALS)

    with pytest.raises(ValueError) as refused:
        calc_kinship(open_vcf(path))

    said = str(refused.value)
    assert "'first' and 'third'" in said
    assert "1 variant is called in the first and 1 in the second" in said
    assert str(path) in said
    assert "position" not in said


def test_the_two_of_the_pair_are_named_as_the_individuals_argument_named_them(
    vcf_of,
) -> None:
    """The same three individuals, with `individuals` naming them in another
    order.

    The core says where the two are among the individuals of the kinship,
    which with `individuals` given is the order of that argument and not the
    order of the file, so a message built from the file's order would name
    the wrong individuals.
    """
    path = vcf_of(NEVER_CALLED_TOGETHER_GTS, NEVER_CALLED_TOGETHER_INDIVIDUALS)

    with pytest.raises(ValueError) as refused:
        calc_kinship(open_vcf(path), individuals=("third", "second", "first"))

    assert "'third' and 'first'" in str(refused.value)


def test_an_individual_with_no_called_genotype_is_named_on_its_own(vcf_of) -> None:
    """One individual of three whose genotype is missing at every variant.

    Its own entry would be divided by no variant at all, so the two of the
    pair are one individual, and the message names that one instead of
    telling the user to leave one of the two out.
    """
    path = vcf_of(
        [
            [[-1, -1], [0, 0], [1, 1]],
            [[-1, -1], [0, 0], [1, 1]],
        ],
        NEVER_CALLED_TOGETHER_INDIVIDUALS,
    )

    with pytest.raises(ValueError) as refused:
        calc_kinship(open_vcf(path))

    said = str(refused.value)
    assert "'first' has no called genotype" in said
    assert "leave it out" in said


def test_a_dataset_in_which_no_variant_varies_is_refused(vcf_of) -> None:
    """Two variants at which every individual holds the same genotype, which
    is what pyNei raises "there is no kinship" for."""
    variants = open_vcf(vcf_of([[[0, 0], [0, 0]], [[1, 1], [1, 1]]]))

    with pytest.raises(ValueError, match="no kinship"):
        calc_kinship(variants)


def test_an_individual_that_is_not_in_the_variants_is_named(vcf_of) -> None:
    """A name of `individuals` that no individual of the source has."""
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))

    with pytest.raises(ValueError, match="i9"):
        calc_kinship(variants, individuals=("i0", "i9"))


def test_the_variants_of_another_type_are_refused_by_name() -> None:
    """The path of the VCF instead of what `open_vcf` gives, which is the
    mistake that is easiest to make."""
    with pytest.raises(TypeError, match="open_vcf"):
        calc_kinship(_panel("panel_called"))


def _kinship_of(values, names, num_vars: int = 10) -> Kinship:
    """A `Kinship` built by hand, as a user brings the matrix plink2 or a
    pedigree gave them."""
    return Kinship(
        matrix=pandas.DataFrame(numpy.array(values), index=names, columns=names),
        num_vars=num_vars,
    )


def test_a_matrix_a_user_built_carries_no_counts_of_a_pass() -> None:
    """A `Kinship` built by hand has `pass_stats` of `None`, since no pass
    over a source gave it, and its individuals are the index of its
    matrix."""
    kinship = _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "b"])

    assert kinship.pass_stats is None
    assert kinship.individuals == ("a", "b")
    assert kinship.num_vars == 10


def test_a_matrix_that_is_not_square_is_refused() -> None:
    """A frame of 2 rows and 3 columns, which is the kinship of no set of
    individuals."""
    with pytest.raises(ValueError, match="2 by 3"):
        Kinship(
            matrix=pandas.DataFrame(
                numpy.zeros((2, 3)), index=["a", "b"], columns=["a", "b", "c"]
            ),
            num_vars=10,
        )


def test_a_matrix_whose_two_sides_are_other_individuals_is_refused() -> None:
    """The rows sorted and the columns not, which leaves every cell the
    kinship of another pair."""
    matrix = pandas.DataFrame(numpy.zeros((2, 2)), index=["a", "b"], columns=["b", "a"])

    with pytest.raises(ValueError, match="the same individuals in the same order"):
        Kinship(matrix=matrix, num_vars=10)


def test_a_matrix_that_is_not_symmetric_names_the_pair_and_its_two_cells() -> None:
    """One cell of a pair filled and the other left at 0, which is what
    building a frame from half a matrix gives."""
    with pytest.raises(ValueError, match="'a'.*'b'"):
        _kinship_of([[1.0, 0.5], [0.0, 1.0]], ["a", "b"])


def test_a_matrix_symmetric_within_the_tolerance_is_taken() -> None:
    """Two cells of a pair that differ by 1e-10 of the largest entry, which
    is inside the 1e-9 of the spec, and two that differ by 1e-8, which is
    outside it."""
    inside = _kinship_of([[1.0, 0.5], [0.5 + 1e-10, 1.0]], ["a", "b"])

    assert inside.individuals == ("a", "b")

    with pytest.raises(ValueError, match="symmetric within"):
        _kinship_of([[1.0, 0.5], [0.5 + 1e-8, 1.0]], ["a", "b"])


def test_filter_individuals_takes_the_rows_and_the_columns_of_those_named() -> None:
    """Three individuals of a kinship of four, in another order than the
    matrix has them.

    It takes entries out and calculates nothing, so `num_vars` and the counts
    of the pass are those of the kinship it came from.
    """
    kinship = calc_kinship(open_vcf(_panel("panel_called")))
    some = ("s004", "s000", "s001")
    filtered = kinship.filter_individuals(some)

    assert filtered.individuals == some
    assert filtered.num_vars == kinship.num_vars
    assert filtered.pass_stats == kinship.pass_stats
    for one in some:
        for other in some:
            assert filtered.matrix.loc[one, other] == kinship.matrix.loc[one, other]


def test_filter_individuals_names_an_individual_that_is_not_in_the_matrix() -> None:
    """A name that the matrix does not hold, which pyNei raises for as
    well."""
    kinship = _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "b"])

    with pytest.raises(ValueError, match="'c'"):
        kinship.filter_individuals(["a", "c"])


def test_a_matrix_that_is_no_frame_names_the_type_it_was_given() -> None:
    """The entries as a numpy array, which has no index to read the names of
    the individuals from."""
    with pytest.raises(TypeError, match="ndarray"):
        Kinship(matrix=numpy.zeros((2, 2)), num_vars=10)


def test_a_matrix_that_holds_what_is_no_number_names_the_argument() -> None:
    """A frame of text, which numpy refuses with words that name neither the
    argument nor the kinship."""
    matrix = pandas.DataFrame(
        [["x", "y"], ["y", "x"]], index=["a", "b"], columns=["a", "b"]
    )

    with pytest.raises(ValueError, match="`matrix` holds what is no number"):
        Kinship(matrix=matrix, num_vars=10)


def test_an_individual_that_is_in_the_matrix_twice_is_refused() -> None:
    """One name on two rows, which makes `filter_individuals` and every
    lookup by that name give two rows where one was asked for."""
    with pytest.raises(ValueError, match="'a' is named twice"):
        _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "a"])


def test_filter_individuals_refuses_a_name_it_was_given_twice() -> None:
    """One name written twice in the call, which would give a 2 x 2 matrix
    for one individual asked for."""
    kinship = _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "b"])

    with pytest.raises(ValueError, match="'a' is named twice"):
        kinship.filter_individuals(["a", "a"])


def test_a_matrix_holding_a_value_that_is_not_finite_names_its_cell() -> None:
    """A NaN in the cell of a pair, which is what a matrix of pyNei holds for
    a pair with no variant called in both of its individuals.

    The matrix is symmetric, NaN in both cells of the pair, so a check that
    took the NaN for an asymmetry would say that the two cells differ by
    `nan` and name no cause a user can act on.
    """
    with pytest.raises(ValueError) as refused:
        _kinship_of(
            [[1.0, numpy.nan], [numpy.nan, 1.0]],
            ["a", "b"],
        )

    said = str(refused.value)
    assert "the cell of the row 'a' and the column 'b' holds nan" in said
    assert "differ by" not in said


def test_a_matrix_holding_an_infinity_is_refused_as_well() -> None:
    """An infinity on the diagonal, which is a number no kinship holds and
    which every value of the matrix would be compared against."""
    with pytest.raises(ValueError, match="holds inf"):
        _kinship_of([[numpy.inf, 0.5], [0.5, 1.0]], ["a", "b"])


def test_filter_individuals_refuses_a_call_that_names_no_individual() -> None:
    """An empty list, which would give a kinship of nobody carrying the
    `num_vars` of the panel it came from."""
    kinship = _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "b"])

    with pytest.raises(ValueError, match="names no individual"):
        kinship.filter_individuals([])


def test_one_name_where_a_sequence_of_them_is_meant_is_refused() -> None:
    """A string, which is a sequence of its letters: `filter_individuals("ab")`
    asks for the individuals `a` and `b`, and for a name of one letter it
    would give a kinship of that individual and say nothing."""
    kinship = _kinship_of([[1.0, 0.5], [0.5, 1.0]], ["a", "b"])

    with pytest.raises(TypeError, match='individuals=\\("a",\\)'):
        kinship.filter_individuals("a")


def test_one_name_where_a_sequence_of_them_is_meant_is_refused_by_calc_kinship(
    vcf_of,
) -> None:
    """The same string given to `calc_kinship`, which without the refusal
    asked the core for the individuals `i`, `0`."""
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))

    with pytest.raises(TypeError, match='individuals=\\("i0",\\)'):
        calc_kinship(variants, individuals="i0")


def test_individuals_that_are_no_sequence_of_names_are_refused(vcf_of) -> None:
    """A whole number where the names are meant, which Python refuses with
    `'int' object is not iterable` and names neither the argument nor the
    call."""
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))

    with pytest.raises(TypeError, match="`individuals` is a sequence"):
        calc_kinship(variants, individuals=4)


def test_calc_kinship_refuses_a_call_that_names_no_individual(vcf_of) -> None:
    """An empty sequence of individuals, which is a kinship of nobody."""
    variants = open_vcf(vcf_of(WORKED_EXAMPLE_GTS, WORKED_EXAMPLE_INDIVIDUALS))

    with pytest.raises(ValueError, match="names no individual"):
        calc_kinship(variants, individuals=())


def test_a_ploidy_above_254_is_refused(vcf_of) -> None:
    """Two individuals whose genotypes hold 255 alleles, which the VCF reader
    takes, its largest ploidy being 255, and which the pass that turns a
    variant into dosages refuses: it writes the genotype of an individual as
    one byte, a dosage from 0 to the ploidy or the code of a genotype with an
    allele missing."""
    genotypes = numpy.zeros((1, 2, 255), dtype=int)
    genotypes[0, 1, :] = 1

    with pytest.raises(ValueError, match="at a ploidy of 255"):
        calc_kinship(open_vcf(vcf_of(genotypes), ploidy=255))


def test_more_individuals_than_the_matrix_counts_in_are_refused(vcf_of) -> None:
    """46341 individuals, one more than the 46340 whose square is inside the
    2147483647 values the linear algebra counts in.

    The check is at the entry of the calculation, before a block is read, so
    the file holds one variant and is 500 kB of text.
    """
    genotypes = numpy.zeros((1, TOO_MANY_INDIVIDUALS, 2), dtype=int)
    genotypes[0, :, 1] = 1

    with pytest.raises(ValueError, match=f"it has {TOO_MANY_INDIVIDUALS} individuals"):
        calc_kinship(open_vcf(vcf_of(genotypes)))
