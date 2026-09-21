"""It writes the positions of the variants that bcftools keeps of each filter.

Run from the root of the repository, with bcftools 1.24 in the PATH, the
version that wrote the files that are committed and the only one the script
runs with:

    uv run --no-project python tests/reference/filters/make_reference.py

It reads tests/reference/vcf/many.vcf, 500 variants of 50 diploid individuals
with 257 half called genotypes, and runs on it the three commands of "How it
is verified" of docs/specs/filters.md, the missing data one, the major allele
frequency one and the observed heterozygosity one, each at the three
thresholds of that table, and the three of them chained at 0.04, 0.8 and 0.5,
in that order. Every variant is given to bcftools, those that failed their
FILTER among them, as popnei reads the file with only_passed false.

It writes, beside itself, one file for each of the eleven sets of filters, with
the position of every variant that set keeps, one position per line, in the
order of the file. The name of a file is `<kind>_<threshold>` for each of its
filters, in the order they were applied, joined by `+`, where the kind is the
one popnei gives the filter in `pass_stats.filtering`:

    missing_data_0.txt   missing_data_0.04.txt   missing_data_0.1.txt
    maf_0.5.txt          maf_0.8.txt             maf_0.95.txt
    obs_het_0.1.txt      obs_het_0.25.txt        obs_het_0.5.txt
    missing_data_0.04+maf_0.8.txt
    missing_data_0.04+maf_0.8+obs_het_0.5.txt

So the number of lines of a file is the number of variants its filters keep:
the nine of the spec's table, and, for the chain, 215 in
missing_data_0.04.txt, which is its first step and one of the nine, then 163
and 106. A cargo test or a pytest test reads the file of the filters it puts
on the variants and compares its positions with the ones popnei gives, and a
test that wants literals copies the first lines of it.

The script checks what it got against the numbers the spec's table gives, the
eleven counts and the first positions of the four rows that have them, and
stops at the first one that differs. The positions of a file are all different,
which is what lets a test compare two sets of variants by them.
"""

import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
MANY_VCF = HERE.parent / "vcf" / "many.vcf"
BCFTOOLS_VERSION = "1.24"

# The kind of each filter and the thresholds of the table of "How it is
# verified" of docs/specs/filters.md, each filter alone on the whole file.
SINGLE_FILTERS = [
    ("missing_data", "0"),
    ("missing_data", "0.04"),
    ("missing_data", "0.1"),
    ("maf", "0.5"),
    ("maf", "0.8"),
    ("maf", "0.95"),
    ("obs_het", "0.1"),
    ("obs_het", "0.25"),
    ("obs_het", "0.5"),
]

# The three filters one after another, each reading what the one before it
# kept, the chain of "How it is verified" of the counts of the same spec.
CHAINED_FILTERS = [("missing_data", "0.04"), ("maf", "0.8"), ("obs_het", "0.5")]

# What the spec says each set of filters keeps of the 500 variants: how many,
# and the first positions where the spec gives them.
EXPECTED = {
    "missing_data_0": (26, [1259, 2110, 2480, 3072, 3257]),
    "missing_data_0.04": (215, []),
    "missing_data_0.1": (455, []),
    "maf_0.5": (35, [1074, 1296, 1481, 1962, 2110]),
    "maf_0.8": (384, []),
    "maf_0.95": (480, []),
    "obs_het_0.1": (22, [1185, 3516, 3923, 4515, 5921]),
    "obs_het_0.25": (79, []),
    "obs_het_0.5": (369, []),
    "missing_data_0.04+maf_0.8": (163, []),
    "missing_data_0.04+maf_0.8+obs_het_0.5": (106, [1111, 1407, 1518]),
}


def check_bcftools():
    """It stops unless the bcftools of the PATH is the version of the files."""
    try:
        version_lines = subprocess.run(
            ["bcftools", "--version"], capture_output=True, text=True, check=True
        ).stdout.splitlines()
    except FileNotFoundError:
        raise SystemExit(
            f"there is no bcftools in the PATH, and this script needs "
            f"bcftools {BCFTOOLS_VERSION}"
        ) from None
    # The first line is "bcftools 1.24".
    version = version_lines[0].removeprefix("bcftools ").strip()
    if version != BCFTOOLS_VERSION:
        raise SystemExit(
            f"the bcftools of the PATH is {version}, and this script needs "
            f"bcftools {BCFTOOLS_VERSION}"
        )


def filter_options(kind, threshold):
    """The options of the `bcftools view` of one filter at one threshold.

    The three commands of "How it is verified" of docs/specs/filters.md, where
    F_MISSING is the fraction of the individuals whose genotype is missing,
    `-Q $t:major` keeps the variants whose major allele frequency is at most
    the threshold, and the third compares the heterozygous individuals with the
    threshold times the called ones, leaving out the variants with nothing
    called.
    """
    if kind == "missing_data":
        return ["-i", f"F_MISSING<={threshold}"]
    if kind == "maf":
        return ["-Q", f"{threshold}:major"]
    if kind == "obs_het":
        few_het = f'N_PASS(GT="het") <= {threshold} * (N_SAMPLES - N_MISSING)'
        something_called = "N_MISSING < N_SAMPLES"
        return ["-i", f"{few_het} && {something_called}"]
    raise ValueError(f"there is no bcftools command for the filter {kind}")


def name_of(filters):
    """The name that the file of a set of filters has, without its suffix."""
    return "+".join(f"{kind}_{threshold}" for kind, threshold in filters)


def view(source_vcf, options, kept_vcf=None):
    """It runs one `bcftools view` and gives the positions of what it keeps.

    With `kept_vcf` it writes what it keeps there, as a VCF with its header,
    for the next filter of a chain to read.
    """
    command = ["bcftools", "view"]
    if kept_vcf is None:
        command.append("-H")
    command += options
    command.append(str(source_vcf))
    kept = subprocess.run(command, capture_output=True, text=True, check=True).stdout
    if kept_vcf is not None:
        kept_vcf.write_text(kept)
    return [
        int(line.split("\t")[1])
        for line in kept.splitlines()
        if not line.startswith("#")
    ]


def run_filters():
    """The positions each set of filters keeps, by the name of its file."""
    kept = {
        name_of([(kind, threshold)]): view(MANY_VCF, filter_options(kind, threshold))
        for kind, threshold in SINGLE_FILTERS
    }
    with tempfile.TemporaryDirectory() as tmp_dir:
        source_vcf = MANY_VCF
        applied = []
        for step, (kind, threshold) in enumerate(CHAINED_FILTERS):
            applied.append((kind, threshold))
            is_last = step == len(CHAINED_FILTERS) - 1
            kept_vcf = None if is_last else Path(tmp_dir) / f"step{step}.vcf"
            positions = view(source_vcf, filter_options(kind, threshold), kept_vcf)
            name = name_of(applied)
            # The first step of the chain is one of the single filters, and
            # the same command run twice has to keep the same variants.
            assert kept.get(name, positions) == positions, name
            kept[name] = positions
            source_vcf = kept_vcf
    return kept


def check(kept):
    """It compares what bcftools kept with the numbers of the spec's table."""
    assert set(kept) == set(EXPECTED), sorted(set(kept) ^ set(EXPECTED))
    for name, (num_vars_kept, first_positions) in EXPECTED.items():
        positions = kept[name]
        assert len(positions) == num_vars_kept, (name, len(positions))
        assert positions[: len(first_positions)] == first_positions, (
            name,
            positions[:5],
        )
        assert len(set(positions)) == len(positions), name


def write(kept):
    for name, positions in kept.items():
        text = "".join(f"{position}\n" for position in positions)
        (HERE / f"{name}.txt").write_text(text)


if __name__ == "__main__":
    check_bcftools()
    kept = run_filters()
    check(kept)
    write(kept)
