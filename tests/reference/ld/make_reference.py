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
    """The variants of `kept` that were kept before `var` in the file, are on
    its chromosome and are no more than `max_dist` base pairs behind it."""
    return [
        before
        for before in kept
        if before < var
        and chroms[before] == chroms[var]
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
    variants or pairs that break it and the number of cases looked at."""
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
        # What a dropped variant is asked about are the variants kept before
        # it, which is what the spec says and what the filter had it to
        # compare with; the `before < var` of `the_window_of` is what leaves
        # the later ones out. Asking about the whole kept set is a weaker
        # check, because a variant kept after the dropped one is behind it by
        # a negative number of base pairs, which is under every max_dist: at
        # max_dist 50000 and max_allowed_r2 0.3, 336 of the 347 dropped
        # variants have a variant kept after them inside that distance, and
        # on the set that the first ten kept variants were taken out of, the
        # weaker check names 45 variants dropped for no reason where this one
        # names 64.
        window = the_window_of(var, kept, max_dist)
        if not any(r2[before, var] > max_allowed_r2 for before in window):
            dropped_with_nothing_above += 1
    return (
        (of_one_dosage, len(kept)),
        (pairs_above, pairs_in_the_window),
        (dropped_with_nothing_above, dropped_with_two_dosages),
    )


def the_lines_of(properties):
    """The three properties written out, a line each."""
    of_one_dosage, above_the_threshold, for_no_reason = properties
    return (
        f"  kept with fewer than two dosages: {of_one_dosage[0]} of "
        f"{of_one_dosage[1]}\n"
        f"  kept pairs inside the window above the threshold: "
        f"{above_the_threshold[0]} of {above_the_threshold[1]}\n"
        f"  dropped with two dosages and nothing above the threshold kept "
        f"before them inside their window: {for_no_reason[0]} of "
        f"{for_no_reason[1]}\n"
    )


def the_sets_written_to(path):
    """The variants that `path` names at each setting, by their index in
    ld.vcf and in the order they are written there."""
    where = {(chroms[var], positions[var]): var for var in range(num_vars)}
    sets = {}
    with open(path) as fh:
        columns = fh.readline().rstrip("\n").split("\t")
        if columns != ["max_dist", "max_allowed_r2", "chrom", "pos"]:
            raise SystemExit(
                f"{path} begins with the columns {columns} and not with "
                "max_dist, max_allowed_r2, chrom and pos"
            )
        for number, line in enumerate(fh, start=2):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 4:
                raise SystemExit(
                    f"line {number} of {path} has {len(fields)} fields and "
                    f"not 4: {line!r}"
                )
            max_dist, max_allowed_r2, chrom, pos = fields
            var = where.get((chrom, int(pos)))
            if var is None:
                raise SystemExit(
                    f"line {number} of {path} names {chrom} {pos}, which is "
                    "not a variant of ld.vcf"
                )
            written = sets.setdefault((int(max_dist), float(max_allowed_r2)), [])
            if written and var <= written[-1]:
                raise SystemExit(
                    f"line {number} of {path} names {chrom} {pos}, which does "
                    "not come after the variant of the line before it in "
                    "ld.vcf, and the properties read these variants in order"
                )
            written.append(var)
    return sets


def against_the_last_kept_variant_alone(max_allowed_r2):
    """The variants kept when a candidate is compared with the last variant
    kept and with nothing else, and neither the chromosome nor the position is
    read, which is the rule of `_filter_chunk_by_ld` of pyNei. The rest is the
    rule above, so a variant of one dosage is dropped and a pair with no r2
    does not drop the candidate."""
    kept, last_kept = [], None
    for var in range(num_vars):
        if not has_two_dosages[var]:
            continue
        if last_kept is None or not (r2[last_kept, var] > max_allowed_r2):
            kept.append(var)
            last_kept = var
    return kept


with open(SP + "/ld.filtered.tsv", "w") as kept_fh:
    kept_fh.write("max_dist\tmax_allowed_r2\tchrom\tpos\n")
    for max_dist, max_allowed_r2 in THE_SETTINGS:
        for var in the_variants_kept(max_dist, max_allowed_r2):
            kept_fh.write(
                f"{max_dist}\t{max_allowed_r2}\t{chroms[var]}\t{positions[var]}\n"
            )

# The properties are worked out over the variants read back from
# ld.filtered.tsv and not over the lists the loop above built, because
# ld.filtered.tsv is the file that the cargo test at
# LdFilteredReader::next_block compares popnei's kept set with: a property
# worked out over the file catches a file that holds something else than the
# rule gives.
the_sets = the_sets_written_to(SP + "/ld.filtered.tsv")
for setting in THE_SETTINGS:
    if setting not in the_sets:
        raise SystemExit(
            f"ld.filtered.tsv names no variant at max_dist {setting[0]} and "
            f"max_allowed_r2 {setting[1]}"
        )

# Each of the three properties is worked out over a set built with the rule
# that property states, so on the four sets above all three are 0 whatever
# that rule is, and running them there cannot fail. Each property is given
# below a set that breaks it, at one setting, so that a run says what the
# check catches. The number named for each of those sets is asserted not to
# be 0, and it is written into ld.filter.properties.txt, which run_plink2.sh
# compares byte for byte with the copy stored in the repository.
DAMAGED_MAX_DIST, DAMAGED_MAX_ALLOWED_R2 = 50000, 0.3
THE_ORDINALS = ("first", "second", "third")


def the_damaged_sets():
    """For each property, by its place in the three, a set that breaks it,
    what that set is and why the number comes out as it does."""
    kept = the_sets[(DAMAGED_MAX_DIST, DAMAGED_MAX_ALLOWED_R2)]
    return (
        (
            0,
            list(range(num_vars)),
            (
                "every variant of ld.vcf kept, the set of a filter that never "
                "asked whether a variant has a second dosage"
            ),
            (
                "the 68 variants whose called genotypes hold one dosage are "
                "in the set, and of the three properties only the first looks "
                "at them"
            ),
        ),
        (
            1,
            against_the_last_kept_variant_alone(DAMAGED_MAX_ALLOWED_R2),
            (
                "each variant compared with the last variant kept and with "
                "nothing else, neither chromosome nor position read, the rule "
                "of _filter_chunk_by_ld of pyNei"
            ),
            (
                "a variant that repeats the variant two before it is kept "
                "whenever the variant in between hid it, so pairs above the "
                "threshold are left inside the window"
            ),
        ),
        (
            2,
            kept[10:],
            (
                f"the {len(kept)} variants kept at these settings with the "
                "first ten of them taken out, the set of a filter that "
                "dropped ten variants it could have kept"
            ),
            (
                "each of the ten, and each variant that only they had "
                "dropped, has two dosages and nothing above the threshold "
                "kept before it inside its window"
            ),
        ),
    )


with open(SP + "/ld.filter.properties.txt", "w") as properties_fh:
    properties_fh.write(
        f"{int((~has_two_dosages).sum())} of the {num_vars} variants of ld.vcf "
        "have one dosage at most among their called genotypes, which the "
        "diagonal of the r2 matrix of plink2 says.\n"
    )
    for max_dist, max_allowed_r2 in THE_SETTINGS:
        kept = the_sets[(max_dist, max_allowed_r2)]
        properties_fh.write(
            f"max_dist {max_dist}, max_allowed_r2 {max_allowed_r2}: "
            f"{len(kept)} variants kept of {num_vars}\n"
            + the_lines_of(the_properties_of(kept, max_dist, max_allowed_r2))
        )
        print(
            f"the filter at {max_dist} bp and {max_allowed_r2} keeps "
            f"{len(kept)} variants of {num_vars}"
        )
    properties_fh.write(
        "\nEach of the four sets above was built with the rule the three "
        "properties state, so there all three are 0 by construction. Below, "
        "each property is given a set that breaks it, at max_dist "
        f"{DAMAGED_MAX_DIST} and max_allowed_r2 {DAMAGED_MAX_ALLOWED_R2}. If "
        "one of these three numbers came out 0, that property would no longer "
        "catch what it is there for, and the program that writes this file "
        "stops rather than write it.\n"
    )
    for which, damaged, what_it_is, why in the_damaged_sets():
        properties = the_properties_of(
            damaged, DAMAGED_MAX_DIST, DAMAGED_MAX_ALLOWED_R2
        )
        broken = properties[which][0]
        if broken == 0:
            raise SystemExit(
                f"the {THE_ORDINALS[which]} property of the item 'The filter "
                "by linkage disequilibrium' of docs/specs/filters.md names 0 "
                f"variants or pairs against {what_it_is}, a set of "
                f"{len(damaged)} variants, at max_dist {DAMAGED_MAX_DIST} and "
                f"max_allowed_r2 {DAMAGED_MAX_ALLOWED_R2}, where it has to "
                f"name some, because {why}"
            )
        properties_fh.write(
            f"The {THE_ORDINALS[which]} property against {what_it_is}: "
            f"{len(damaged)} of the {num_vars} variants kept. The "
            f"{THE_ORDINALS[which]} number below is {broken} because {why}.\n"
            + the_lines_of(properties)
        )
        print(
            f"the {THE_ORDINALS[which]} property names {broken} on the set "
            f"of {len(damaged)} variants that breaks it"
        )
