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

# The variants that the filter by linkage disequilibrium of
# docs/specs/filters.md keeps of ld.vcf, at the four settings of the table of
# "How it is verified" of the item "The filter by linkage disequilibrium",
# and the three properties that item asks of the set kept.
#
# The rule is stated here a second time, in Python and over the r2 that
# plink2 wrote into ld.unphased.vcor2.bin, so that the set does not depend on
# popnei's own arithmetic: a variant is kept when its called genotypes hold
# two dosages at least and its r2 against every variant kept no more than
# max_dist base pairs behind it on its chromosome is at most max_allowed_r2,
# and a pair with no r2 does not drop it. What ties popnei to this set is a
# cargo test at LdFilteredReader::next_block, which asserts that the reader
# gives these chromosomes and positions and no others.
#
# The diagonal of plink2's matrix is the r2 of a variant with itself, 1 for a
# variant whose called genotypes hold two dosages at least and NaN for one
# whose genotypes all hold the same dosage, so it is where a variant of one
# dosage is read from: of the 500 variants, 68 have no second dosage.
THE_SETTINGS = ((10000, 0.1), (10000, 0.3), (50000, 0.3), (250000, 0.3))

chroms, positions, identifiers = [], [], []
with open(SP + "/ld.vcf") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.split("\t", 3)
        chroms.append(fields[0])
        positions.append(int(fields[1]))
        identifiers.append(fields[2])

# The matrix read is the one stored in the repository, because plink2 runs
# after this script and has not written the new one yet; run_plink2.sh is
# what compares the stored matrix with the one plink2 writes, and the rows
# of the matrix are named in the file beside it, in their order.
stored_matrix = os.path.join(HERE, "ld.unphased.vcor2.bin")
with open(stored_matrix + ".vars") as fh:
    rows_of_the_matrix = fh.read().split()
if rows_of_the_matrix != identifiers:
    raise SystemExit(
        f"ld.vcf has {len(identifiers)} variants and {stored_matrix}.vars names "
        f"{len(rows_of_the_matrix)}, and they are not the same variants in the "
        "same order, so the matrix is not the one of this VCF"
    )
num_vars = len(identifiers)
r2 = numpy.fromfile(stored_matrix, dtype="<f8").reshape(num_vars, num_vars)
has_two_dosages = ~numpy.isnan(numpy.diagonal(r2))


def the_window_of(var, kept, max_dist):
    """The variants of `kept` that are on the chromosome of `var` and no
    more than `max_dist` base pairs behind it."""
    return [
        before
        for before in kept
        if chroms[before] == chroms[var]
        and positions[var] - positions[before] <= max_dist
    ]


def the_variants_kept(max_dist, max_allowed_r2):
    """The variants kept, by their index in the file, in the order they come."""
    kept = []
    for var in range(num_vars):
        if not has_two_dosages[var]:
            continue
        window = the_window_of(var, kept, max_dist)
        # An r2 that is not a number is above no threshold, so `not (r2 >
        # threshold)` keeps the candidate where `r2 <= threshold` would drop
        # it.
        if all(not (r2[before, var] > max_allowed_r2) for before in window):
            kept.append(var)
    return kept


def the_properties_of(kept, max_dist, max_allowed_r2):
    """The three properties of "How it is verified", each as the number of
    variants or pairs that break it, with the number of cases looked at."""
    of_one_dosage = sum(1 for var in kept if not has_two_dosages[var])
    pairs_in_the_window, pairs_above = 0, 0
    for at, var in enumerate(kept):
        for before in the_window_of(var, kept[:at], max_dist):
            pairs_in_the_window += 1
            if r2[before, var] > max_allowed_r2:
                pairs_above += 1
    dropped_with_two_dosages, dropped_with_nothing_above = 0, 0
    was_kept = set(kept)
    for var in range(num_vars):
        if var in was_kept or not has_two_dosages[var]:
            continue
        dropped_with_two_dosages += 1
        window = the_window_of(var, kept, max_dist)
        if not any(r2[before, var] > max_allowed_r2 for before in window):
            dropped_with_nothing_above += 1
    return (
        f"max_dist {max_dist}, max_allowed_r2 {max_allowed_r2}: "
        f"{len(kept)} variants kept of {num_vars}\n"
        f"  kept with fewer than two dosages: {of_one_dosage} of {len(kept)}\n"
        f"  kept pairs inside the window above the threshold: "
        f"{pairs_above} of {pairs_in_the_window}\n"
        f"  dropped with two dosages and nothing above the threshold kept "
        f"before them inside their window: {dropped_with_nothing_above} of "
        f"{dropped_with_two_dosages}\n"
    )


with (
    open(SP + "/ld.filtered.tsv", "w") as kept_fh,
    open(SP + "/ld.filter.properties.txt", "w") as properties_fh,
):
    kept_fh.write("max_dist\tmax_allowed_r2\tchrom\tpos\n")
    properties_fh.write(
        f"{int((~has_two_dosages).sum())} of the {num_vars} variants of ld.vcf "
        "have one dosage at most among their called genotypes, which the "
        "diagonal of the r2 matrix of plink2 says.\n"
    )
    for max_dist, max_allowed_r2 in THE_SETTINGS:
        kept = the_variants_kept(max_dist, max_allowed_r2)
        for var in kept:
            kept_fh.write(
                f"{max_dist}\t{max_allowed_r2}\t{chroms[var]}\t{positions[var]}\n"
            )
        properties_fh.write(the_properties_of(kept, max_dist, max_allowed_r2))
        print(
            f"the filter at {max_dist} bp and {max_allowed_r2} keeps "
            f"{len(kept)} variants of {num_vars}"
        )
