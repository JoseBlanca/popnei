"""It writes a VCF of 1000 individuals and as many variants as it is asked
for.

Two benchmarks read what it writes: the VCF of 100000 variants that task 5.2
of docs/plans/vcf-to-blocks.md times the VCF reader on, which is what it
writes when no number of variants is given, and the VCF of 20000 variants,
the panel that "The compression" of docs/specs/io_vars.md measured, which is
written as a vars file and given to `vars_file.rs`.

The genotypes are simulated by `simulate_genotypes` of
test/gwas_reference/make_reference.py of pyNei, copied here with the number
of variants as an argument in the place of its constant, with NUM_SAMPLES
1000 in the place of the 200 of that file, which is the panel section 6 of
docs/rust_core.md names. The VCF is written by `write_vcf` of the same file:
the same header, the same nine first columns with `.` in QUAL, FILTER and
INFO, `A` and `T` as the alleles, and `./.` where a genotype is missing. Here
it writes the same bytes through numpy instead of one f-string per genotype,
because 100 million f-strings take minutes.

What differs from that script: it stops after the genotypes. It does not
simulate the traits, write the vars file or run plink2 and R, so the missing
mask is drawn from the generator right after the genotypes instead of after
the traits, and the missing genotypes are not the ones a full run would have
drawn. Everything else, the seed included, is the script's.

Every draw of the generator is of the number of variants it is given, so a
VCF of one number of variants is not a prefix of a VCF of another: the two
files hold different genotypes.

    uv run --no-project --with numpy python make_big_vcf.py <out.vcf> [num_vars]
"""

import sys

import numpy

SEED = 42
NUM_SAMPLES = 1000
# How many variants it writes when the command line does not say.
DEFAULT_NUM_VARS = 100000
NUM_CHROMS = 2
NUM_POPS = 3
FAMILY_SIZE = 4
FST = 0.1
MISSING_RATE = 0.03

VARS_PER_WRITE = 1000


def simulate_genotypes(rng, num_vars):
    num_families = NUM_SAMPLES // FAMILY_SIZE
    family_pops = rng.integers(0, NUM_POPS, size=num_families)
    p_anc = rng.uniform(0.1, 0.9, size=num_vars)
    a = p_anc * (1 - FST) / FST
    b = (1 - p_anc) * (1 - FST) / FST
    p_pop = numpy.stack([rng.beta(a, b) for _ in range(NUM_POPS)])
    alleles = numpy.empty((num_vars, NUM_SAMPLES, 2), dtype=numpy.int8)
    pops = numpy.repeat(family_pops, FAMILY_SIZE)
    for family_idx, pop in enumerate(family_pops):
        parents = (
            rng.uniform(size=(2, num_vars, 2)) < p_pop[pop][None, :, None]
        ).astype(numpy.int8)
        for child_idx in range(FAMILY_SIZE):
            sample_idx = family_idx * FAMILY_SIZE + child_idx
            for parent_idx in range(2):
                picked = rng.integers(0, 2, size=num_vars)
                alleles[:, sample_idx, parent_idx] = parents[
                    parent_idx, numpy.arange(num_vars), picked
                ]
    return alleles, pops


# The five texts a genotype of a biallelic variant takes, at the place of
# the code `2 * first_allele + second_allele`, and the missing one last.
GT_TEXTS = numpy.frombuffer(b"0/00/11/01/1./.", dtype=numpy.uint8).reshape(5, 3)


def write_vcf(path, alleles, samples, chroms, poss, ids):
    with open(path, "wb") as fhand:
        fhand.write(b"##fileformat=VCFv4.2\n")
        for chrom in dict.fromkeys(chroms):
            fhand.write(f"##contig=<ID={chrom}>\n".encode())
        fhand.write(
            b'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        )
        fhand.write(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t".encode()
            + "\t".join(samples).encode()
            + b"\n"
        )
        for start in range(0, alleles.shape[0], VARS_PER_WRITE):
            block = alleles[start : start + VARS_PER_WRITE].astype(numpy.int16)
            codes = 2 * block[:, :, 0] + block[:, :, 1]
            codes[block[:, :, 0] < 0] = 4
            texts = GT_TEXTS[codes]
            # every genotype followed by a tab, and the last one of the line
            # by the end of the line
            line = numpy.empty((texts.shape[0], texts.shape[1], 4), numpy.uint8)
            line[:, :, :3] = texts
            line[:, :, 3] = ord("\t")
            line[:, -1, 3] = ord("\n")
            rows = line.reshape(texts.shape[0], -1)
            for offset in range(rows.shape[0]):
                var_idx = start + offset
                fhand.write(
                    f"{chroms[var_idx]}\t{poss[var_idx]}\t{ids[var_idx]}"
                    "\tA\tT\t.\t.\t.\tGT\t".encode()
                )
                fhand.write(rows[offset].tobytes())


def main():
    out = sys.argv[1]
    num_vars = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_NUM_VARS
    # The variants are split evenly between the chromosomes, so a number
    # that does not divide by them would leave the last ones with no
    # chromosome and no position.
    if num_vars <= 0 or num_vars % NUM_CHROMS:
        raise SystemExit(
            f"{num_vars} variants: the number has to be 1 or more and a "
            f"multiple of the {NUM_CHROMS} chromosomes"
        )
    rng = numpy.random.default_rng(SEED)
    alleles, _pops = simulate_genotypes(rng, num_vars)
    samples = [f"s{idx:03d}" for idx in range(NUM_SAMPLES)]
    vars_per_chrom = num_vars // NUM_CHROMS
    chroms = [
        f"chr{idx + 1}" for idx in range(NUM_CHROMS) for _ in range(vars_per_chrom)
    ]
    poss = [
        1000 * (idx + 1) for _ in range(NUM_CHROMS) for idx in range(vars_per_chrom)
    ]
    ids = [f"var{idx:04d}" for idx in range(num_vars)]
    is_missing = rng.uniform(size=(num_vars, NUM_SAMPLES)) < MISSING_RATE
    alleles[is_missing] = -1
    write_vcf(out, alleles, samples, chroms, poss, ids)
    print(f"{out}: {num_vars} variants x {NUM_SAMPLES} individuals")


if __name__ == "__main__":
    sys.exit(main())
