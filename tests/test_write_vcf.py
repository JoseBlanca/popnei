"""`write_vcf`: the variants of a pass written as a VCF.

The cases are the ones `docs/specs/io_vcf.md` gives to pytest under "How it
is verified" of the writer: `write.vcf` read with every line and written
back, plain and bgzipped, and the file written from the vars file of
`write.vcf`, whose five lines are literals of that spec; a path that a file
is already at, and a VCF with a wrong line after 250 good ones, whose file
is taken away. The filter of individuals that keeps `c` and `a` is
compared with what bcftools 1.24 wrote, `write.c_a.bcftools.vcf`, which
`tests/reference/vcf/make_reference.py` stores.

pyNei is not run here: it has no VCF writer.
"""

import errno
import gzip
from pathlib import Path

import pytest
from popnei import (
    VcfWritten,
    _core,
    open_vars,
    open_vcf,
    write_vars,
    write_vcf,
)

# The line bcftools adds after `##fileformat` and popnei does not.
BCFTOOLS_PASS_LINE = '##FILTER=<ID=PASS,Description="All filters passed">\n'

# The file written from the vars file of `write.vcf` read with the default,
# the header and the five lines of "How it is verified" of the writer. The
# vars file of a VCF keeps whether each variant passed, so the header has
# the line of `FAIL`, which no line here has.
FROM_THE_VARS_FILE = (
    "##fileformat=VCFv4.3\n"
    '##FILTER=<ID=FAIL,Description="It failed a filter of the VCF the variants were read from">\n'
    "##contig=<ID=chr1,length=2000>\n"
    "##contig=<ID=chr2,length=1500>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
    "chr1\t100\trs1\tA\tT\t29.5\t.\t.\tGT\t0/1\t0/1\t1/1\n"
    "chr1\t1000\trs3\tG\tC,T\t50\t.\t.\tGT\t1/2\t0/1\t0/0\n"
    "chr1\t1001\t.\tC\t.\t12\t.\t.\tGT\t0/0\t0/0\t0/0\n"
    "chr2\t1\t.\tT\tG\t40\t.\t.\tGT\t1/1\t0/.\t0/0\n"
    "chr2\t1500\trs6\tA\tG\t33\t.\t.\tGT\t0/0\t1/0\t0/0\n"
)


def test_write_vcf_writes_write_vcf_back_as_it_was_plain_and_bgzipped(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Read with every line and written with no step, the bytes are the
    file's; bgzipped, they are what the members decompress to, and the file
    ends with the empty member of 28 bytes that bgzip writes."""
    source = reference_vcf_dir / "write.vcf"
    variants = open_vcf(source, only_passed=False)

    plain = tmp_path / "written.vcf"
    written = write_vcf(variants, plain)
    assert isinstance(written, VcfWritten)
    assert written.pass_stats.num_vars == 6
    assert plain.read_bytes() == source.read_bytes()

    bgzipped = tmp_path / "written.vcf.gz"
    assert write_vcf(variants, bgzipped).pass_stats.num_vars == 6
    assert gzip.decompress(bgzipped.read_bytes()) == source.read_bytes()
    assert bgzipped.read_bytes()[-28:] == bytes.fromhex(
        "1f8b08040000000000ff0600424302001b0003000000000000000000"
    )


def test_write_vcf_bgzips_a_path_that_ends_in_gz_in_any_case(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """`.GZ` and `.Gz` are read as `.gz`, and `.gzip` is not."""
    source = reference_vcf_dir / "write.vcf"
    variants = open_vcf(source, only_passed=False)
    for name in ("a.VCF.GZ", "b.vcf.Gz"):
        path = tmp_path / name
        write_vcf(variants, path)
        assert gzip.decompress(path.read_bytes()) == source.read_bytes()
    plain = tmp_path / "c.vcf.gzip"
    write_vcf(variants, plain)
    assert plain.read_bytes() == source.read_bytes()


def test_write_vcf_writes_the_five_lines_of_the_vars_file_of_write_vcf(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """Read with the default, written as a vars file and written back.

    The phase of `0|1` and of `1|1` is gone, since a vars file keeps none,
    and FILTER and INFO are a dot."""
    vars_path = tmp_path / "write.vars"
    write_vars(open_vcf(reference_vcf_dir / "write.vcf"), vars_path)
    path = tmp_path / "from_the_vars_file.vcf"

    written = write_vcf(open_vars(vars_path), path)

    assert written.pass_stats.num_vars == 5
    assert path.read_text() == FROM_THE_VARS_FILE


def test_write_vcf_with_the_filter_of_c_and_a_is_what_bcftools_writes(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """AC and AN are taken out of the lines and of the header, as
    `bcftools annotate -x INFO/AC,INFO/AN` takes them out, and the counts of
    the filter are those of a filter that kept every variant."""
    variants = open_vcf(reference_vcf_dir / "write.vcf", only_passed=False)
    variants.filter_individuals(["c", "a"])
    path = tmp_path / "c_a.vcf"

    written = write_vcf(variants, path)

    bcftools = (reference_vcf_dir / "write.c_a.bcftools.vcf").read_text()
    assert BCFTOOLS_PASS_LINE in bcftools
    assert path.read_text() == bcftools.replace(BCFTOOLS_PASS_LINE, "", 1)
    assert written.pass_stats.num_vars == 6


def test_write_vcf_refuses_a_path_that_a_file_is_already_at(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """The second call writes nothing and leaves the first file as it was."""
    variants = open_vcf(reference_vcf_dir / "write.vcf")
    path = tmp_path / "write.vcf.gz"
    write_vcf(variants, path)
    written = path.read_bytes()

    with pytest.raises(ValueError, match="already") as refusal:
        write_vcf(variants, path)

    assert str(refusal.value).startswith(str(path))
    assert path.read_bytes() == written


# How many individuals the VCF of the test of a wrong line has: the writer
# reads a VCF of 5000 individuals in blocks of 200 variants, a fifth of 5
# million genotypes over them, so the 200 lines of its first block, 4 MB of
# text, reach the file before the block of the wrong line is read.
INDIVIDUALS_OF_THE_WRONG_LINE = 5000


def test_write_vcf_leaves_no_file_when_the_vcf_has_a_wrong_line_after_250_good_ones(
    tmp_path: Path,
) -> None:
    """A line of the columns of 4999 individuals, the 251st of a VCF of 5000.

    The first block of 200 lines is written before the error, which is the
    one of the source and names the VCF, and the path the user wrote to is
    free afterwards, plain and bgzipped."""
    individuals = [f"i{index}" for index in range(INDIVIDUALS_OF_THE_WRONG_LINE)]
    genotypes = "\t".join(["0/1"] * INDIVIDUALS_OF_THE_WRONG_LINE)
    lines = [
        "##fileformat=VCFv4.3",
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            + individuals
        ),
    ]
    lines += [
        f"chr1\t{pos}\t.\tA\tT\t.\tPASS\t.\tGT\t{genotypes}" for pos in range(1, 251)
    ]
    # The column of the last individual is missing.
    lines.append("chr1\t251\t.\tA\tT\t.\tPASS\t.\tGT\t" + genotypes[4:])
    vcf_path = tmp_path / "wrong_after_250.vcf"
    vcf_path.write_text("\n".join(lines) + "\n")
    for name in ("half_way.vcf", "half_way.vcf.gz"):
        path = tmp_path / name
        with pytest.raises(ValueError, match="4999 individuals") as refusal:
            write_vcf(open_vcf(vcf_path), path)
        assert str(refusal.value).startswith(str(vcf_path))
        # The header is the lines 1 and 2, so the 251st variant is line 253.
        assert "line 253 " in str(refusal.value)
        assert not path.exists()


def test_write_vcf_gives_the_error_of_the_file_system_for_a_path_of_no_file(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A directory where the VCF goes, and one that is not there, are the
    `IsADirectoryError` and the `FileNotFoundError` of those paths."""
    variants = open_vcf(reference_vcf_dir / "write.vcf")

    with pytest.raises(IsADirectoryError) as refusal:
        write_vcf(variants, tmp_path)
    assert refusal.value.errno == errno.EISDIR
    assert refusal.value.filename == str(tmp_path)

    of_no_directory = tmp_path / "no_such_directory" / "write.vcf.gz"
    with pytest.raises(FileNotFoundError) as refusal:
        write_vcf(variants, of_no_directory)
    assert refusal.value.filename == str(of_no_directory)


def test_write_vcf_says_what_it_takes_when_it_is_given_a_path(
    reference_vcf_dir: Path, tmp_path: Path
) -> None:
    """A user who gives the VCF where its variants go gets a `TypeError`
    that says where the variants come from, and no file is made."""
    path = tmp_path / "write.vcf"

    with pytest.raises(TypeError, match="open_vcf"):
        write_vcf(str(reference_vcf_dir / "write.vcf"), path)  # type: ignore[arg-type]

    assert not path.exists()


def test_what_a_user_reads_of_write_vcf_is_written_in_the_package() -> None:
    """The private module explains nothing; the package is the API."""
    assert _core.write_vcf.__doc__ is None
    assert write_vcf.__doc__ is not None
