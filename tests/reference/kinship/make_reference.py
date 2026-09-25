r"""It writes the datasets of the genomic relationship matrix and what plink2 gives for them.

Run from the root of the repository:

    uv run python tests/reference/kinship/make_reference.py

It needs two things outside the repository, and it says which one is missing
when it cannot find it:

- plink2 v2.0.0-a.7.7, the build the numbers of `docs/specs/kinship.md` were
  taken with, on the PATH. The script refuses any other version, because
  `--make-rel` is what every literal of that spec comes from.
- /Users/jose/devel/pynei/test/gwas_reference/sim.vars, the vars file of pyNei
  that holds the panel with every genotype called. popnei does not read a vars
  file of pyNei, so the script reads it with pyNei, which is a development
  dependency of popnei, and writes it as a VCF.

The two datasets are the ones of "How it is verified" of
`docs/specs/kinship.md`, and they are the same 200 individuals and 1200
biallelic diploid variants twice, once with every genotype called and once
with 3 in 100 of them missing whole:

- `panel_called`, from `sim.vars`, written here as `panel_called.vcf.gz`.
- `panel`, which is `tests/reference/dists/panel.vcf.gz` of
  `docs/specs/dists.md`, already in the repository and not written again. The
  kinship needs it because the denominator of a pair is the variants called in
  both of its individuals, and with nothing missing that denominator is the
  same for every pair and the rule is never exercised.

It writes, beside itself:

- `panel_called.vcf.gz`, the genotypes as a gzipped VCF that popnei's
  `open_vcf` reads. The timestamp of gzip is fixed to 0, so that running the
  script again gives the same bytes.
- `<name>.plink2.rel.gz`, what `plink2 --make-rel square` writes: 200 lines of
  200 numbers separated by tabs, the matrix row after row, gzipped because the
  text is 1 MB. It is what the table of literals of the spec is read from, by
  eye, and it holds six significant digits.
- `<name>.plink2.rel.bin.gz`, what `plink2 --make-rel square bin` writes: the
  same matrix as 40000 little endian `f64`, row after row, gzipped. The cargo
  tests compare against this one, because the text rounds an entry near 1 by
  up to 5e-6 and a comparison with it leaves no room to find an error of the
  arithmetic in.
- `<name>.plink2.rel.id`, the individuals of that matrix in its order, so that
  a test does not have to trust that plink2 kept the order of the VCF.

At the end it checks the literals that `docs/specs/kinship.md` writes into the
cargo tests against what it has just produced, and says so for each one, so
that a run of this script is also a check that the spec has not drifted. It
also says how far the text is from the binary, which is plink2's rounding and
nothing of popnei.
"""

from pathlib import Path
import gzip
import shutil
import subprocess
import sys
import tempfile

import numpy
import pandas

REF_DIR = Path(__file__).parent
DISTS_PANEL = REF_DIR.parent / "dists" / "panel.vcf.gz"
PYNEI_SIM = Path("/Users/jose/devel/pynei/test/gwas_reference/sim.vars")

PLINK2_VERSION = "v2.0.0-a.7.7"

# The literals of "How it is verified" of docs/specs/kinship.md: the dataset,
# the two individuals of the pair, counting from 0, and what plink2 wrote.
LITERALS = [
    ("panel_called", 0, 0, 1.09309),
    ("panel_called", 1, 1, 1.22825),
    ("panel_called", 0, 1, 0.648081),
    ("panel_called", 0, 2, 0.615611),
    ("panel_called", 0, 4, -0.0945533),
    ("panel_called", 0, 199, -0.0760273),
    ("panel_called", 100, 101, 0.604995),
    ("panel", 0, 0, 1.09626),
    ("panel", 0, 1, 0.650379),
    ("panel", 0, 4, -0.103505),
    ("panel", 100, 101, 0.604119),
]


def refuse_a_missing_tool():
    if shutil.which("plink2") is None:
        sys.exit("plink2 is not on the PATH; this script needs " + PLINK2_VERSION)
    printed = subprocess.run(
        ["plink2", "--version"], capture_output=True, text=True, check=True
    ).stdout
    if PLINK2_VERSION not in printed:
        sys.exit(f"this script needs plink2 {PLINK2_VERSION}, and found: {printed.strip()}")
    if not PYNEI_SIM.exists():
        sys.exit(
            f"{PYNEI_SIM} is not there. It is the panel with every genotype called, "
            "and it comes from a checkout of pyNei beside popnei."
        )
    if not DISTS_PANEL.exists():
        sys.exit(f"{DISTS_PANEL} is not there; it is in the repository.")


def vcf_text(alleles, individuals, chroms, poss, ids):
    """The genotypes as VCF text, a missing genotype missing whole.

    Which nucleotide each allele gets changes no entry of the matrix: the
    kinship is built from the dosages, and REF and ALT only name them.
    """
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={chrom}>" for chrom in dict.fromkeys(chroms)]
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(individuals)
    )
    for idx in range(alleles.shape[0]):
        gts = "\t".join(
            "./." if first < 0 else f"{first}/{second}"
            for first, second in alleles[idx].tolist()
        )
        lines.append(
            f"{chroms[idx]}\t{poss[idx]}\t{ids[idx]}\tA\tT\t.\t.\t.\tGT\t{gts}"
        )
    return "\n".join(lines) + "\n"


def write_panel_called():
    """It reads pyNei's vars file of the panel and writes it as a gzipped VCF."""
    from pynei import load_vars

    variants = load_vars(PYNEI_SIM)
    chunks = list(variants.iter_vars_chunks())
    alleles = numpy.concatenate([chunk.gts.gt_values for chunk in chunks])
    info = pandas.concat([chunk.vars_info for chunk in chunks], ignore_index=True)
    text = vcf_text(
        alleles,
        list(variants.samples),
        list(info["chrom"]),
        list(info["pos"]),
        list(info["id"]),
    )
    path = REF_DIR / "panel_called.vcf.gz"
    with (
        path.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed,
    ):
        compressed.write(text.encode())
    return path


def gzipped_copy_of(source, target):
    """`source` gzipped into `target`, with the timestamp of gzip fixed to 0."""
    with (
        source.open("rb") as src,
        target.open("wb") as raw,
        gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as dst,
    ):
        shutil.copyfileobj(src, dst)


def run_plink2(vcf_gz, name):
    """`plink2 --make-rel square` on a gzipped VCF, twice: the text and the bits.

    plink2 reads a gzipped VCF, and writes the matrix as text with six
    significant digits, or as `f64` with `bin`, and the individuals in a file
    of its own. Both runs are made, into `<name>.plink2.rel.gz` and
    `<name>.plink2.rel.bin.gz`: the text is what a reader checks the table of
    the spec against and the bits are what the tests compare within 1e-12
    relative, which the text has no room for.
    """
    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        for modifiers, written, kept in (
            (["square"], f"{name}.rel", f"{name}.plink2.rel.gz"),
            (["square", "bin"], f"{name}.rel.bin", f"{name}.plink2.rel.bin.gz"),
        ):
            subprocess.run(
                ["plink2", "--vcf", str(vcf_gz.resolve()), "--make-rel", *modifiers,
                 "--out", name],
                cwd=work, check=True, capture_output=True,
            )
            gzipped_copy_of(work / written, REF_DIR / kept)
        shutil.copy(work / f"{name}.rel.id", REF_DIR / f"{name}.plink2.rel.id")


def read_matrix(name):
    with gzip.open(REF_DIR / f"{name}.plink2.rel.gz", "rt") as fhand:
        return numpy.loadtxt(fhand)


def read_binary_matrix(name):
    with gzip.open(REF_DIR / f"{name}.plink2.rel.bin.gz", "rb") as fhand:
        values = numpy.frombuffer(fhand.read(), dtype="<f8")
    side = round(values.size**0.5)
    return values.reshape(side, side)


def main():
    refuse_a_missing_tool()
    called = write_panel_called()
    print(f"wrote {called.name}")
    for name, vcf in (("panel_called", called), ("panel", DISTS_PANEL)):
        run_plink2(vcf, name)
        print(
            f"wrote {name}.plink2.rel.gz, {name}.plink2.rel.bin.gz and "
            f"{name}.plink2.rel.id"
        )

    matrices = {name: read_matrix(name) for name in ("panel_called", "panel")}
    for name, matrix in matrices.items():
        individuals = pandas.read_csv(REF_DIR / f"{name}.plink2.rel.id", sep="\t")
        order_kept = list(individuals.iloc[:, -1]) == [
            f"s{idx:03d}" for idx in range(matrix.shape[0])
        ]
        print(f"{name}: {matrix.shape[0]} individuals, plink2 kept the order of "
              f"the VCF: {order_kept}")
        of_the_bits = read_binary_matrix(name)
        if of_the_bits.shape != matrix.shape:
            sys.exit(f"{name}: the text is {matrix.shape} and the bits {of_the_bits.shape}")
        print(f"{name}: the text is at most {numpy.abs(matrix - of_the_bits).max():.3g} "
              "from the bits, which is the rounding of its six digits")

    wrong = 0
    for name, row, col, expected in LITERALS:
        got = matrices[name][row, col]
        ok = abs(got - expected) < 5e-7
        wrong += not ok
        print(f"  {'ok  ' if ok else 'WRONG'} {name}[{row}, {col}] = {got} , "
              f"the spec says {expected}")
    if wrong:
        sys.exit(f"{wrong} of the {len(LITERALS)} literals of the spec do not match")
    print(f"the {len(LITERALS)} literals of docs/specs/kinship.md match")


if __name__ == "__main__":
    sys.exit(main())
