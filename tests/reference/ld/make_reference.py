import os
"""Writes the LD reference VCF: two chromosomes whose r2 decays with distance."""
import numpy
import pynei
SP = os.environ.get("LD_WORK", ".")
NUM_IND, NUM_VAR_PER_CHROM, SPACING = 100, 250, 1000
RECOMB = 0.02          # per 1000 bp between one variant and the next
MISSING_RATE = 0.03
rng = numpy.random.default_rng(29)

def haplotypes(num_haps, num_vars):
    """A pool of 4 founder haplotypes recombined along the chromosome."""
    founders = rng.integers(0, 2, size=(4, num_vars))
    out = numpy.empty((num_haps, num_vars), dtype=numpy.int8)
    for h in range(num_haps):
        which = rng.integers(0, 4)
        for v in range(num_vars):
            if v and rng.random() < RECOMB:
                which = rng.integers(0, 4)
            out[h, v] = founders[which, v]
    return out

chrom_gts = []
for _ in range(2):
    haps = haplotypes(NUM_IND * 2, NUM_VAR_PER_CHROM)
    gts = haps.reshape(NUM_IND, 2, NUM_VAR_PER_CHROM).transpose(2, 0, 1)
    chrom_gts.append(gts)

names = [f"i{i:03d}" for i in range(NUM_IND)]
with open(SP + "/ld.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.4\n##contig=<ID=chr1>\n##contig=<ID=chr2>\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
    idx = 0
    for chrom, gts in zip(("chr1", "chr2"), chrom_gts):
        missing = rng.random((NUM_VAR_PER_CHROM, NUM_IND)) < MISSING_RATE
        for v in range(NUM_VAR_PER_CHROM):
            cells = []
            for i in range(NUM_IND):
                cells.append("./." if missing[v, i] else f"{gts[v, i, 0]}/{gts[v, i, 1]}")
            fh.write(f"{chrom}\t{(v + 1) * SPACING}\tv{idx:04d}\tA\tC\t.\tPASS\t.\tGT\t"
                     + "\t".join(cells) + "\n")
            idx += 1
print("written", idx, "variants")

# Everything above asks numpy.random.default_rng(29) for its numbers, and
# every literal of docs/specs/ld.md and of the filter item of
# docs/specs/filters.md depends on the order in which it asked, so what is
# added to this script goes below this line, where it cannot move them: it
# reads files that are already in the repository and writes numbers beside
# ld.vcf.
HERE = os.path.dirname(os.path.abspath(__file__))
MANY_VCF = os.path.join(HERE, "..", "vcf", "many.vcf")

# The dosages of tests/reference/vcf/many.vcf, 500 variants of 50 diploid
# individuals with 54 variants of more than two alleles and 257 of its 25000
# genotypes half called, as pyNei's to_012 gives them: for each genotype, how
# many of its alleles are not the major allele of its variant, with the
# called allele of a half called genotype counted among the alleles the
# major one is chosen from, and -1 for a genotype that has an allele
# missing and so no dosage. It is the rule of docs/specs/pca.md, which the
# LdDosages of docs/specs/ld.md follows, and the cargo test at
# LdDosages::dosages compares every one of the 25000 with this file.
#
# The whole file is read as one chunk, and the counts of the alleles of a
# variant run over every individual either way: to_012 chooses the major
# allele of each variant on its own.
many = pynei.vars_from_vcf(MANY_VCF, desired_num_vars_per_chunk=1000)
dosages = numpy.vstack([chunk.gts.to_012() for chunk in many.iter_vars_chunks()])
numpy.savetxt(SP + "/many.pynei.dosages.tsv", dosages, fmt="%d", delimiter="\t")
print("written the dosages of", dosages.shape[0], "variants of", dosages.shape[1], "individuals")
