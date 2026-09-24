"""The five per population quantities, computed as the spec defines them."""

import gzip
from math import comb

import numpy as np


def read_vcf(path):
    op = gzip.open if str(path).endswith(".gz") else open
    names, rows = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                names = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            gts = []
            for cell in f[9:]:
                g = cell.split(":")[0].replace("|", "/").split("/")
                gts.append([-1 if a == "." else int(a) for a in g])
            rows.append(gts)
    return names, np.array(rows, dtype=np.int16)


def counts_of(gts_pop, num_alleles):
    """The allele counts of one population at one variant, and the called ones."""
    alleles = gts_pop[gts_pop >= 0]
    counts = np.bincount(alleles, minlength=num_alleles).astype(np.int64)
    return counts, int(counts.sum())


def expected_alleles(counts, c, g):
    """The alleles expected in a draw of `g` of the `c` called, no replacement."""
    total = comb(c, g)
    return sum(1.0 - comb(c - int(n), g) / total for n in counts if n > 0)


def chance_variable(counts, c, g):
    """The chance that a draw of `g` of the `c` called shows more than one allele."""
    total = comb(c, g)
    return 1.0 - sum(comb(int(n), g) / total for n in counts if n > 0)


def chance_in_draw(n, c, g):
    """The chance that an allele called `n` times appears in a draw of `g`."""
    return 1.0 - comb(c - int(n), g) / comb(c, g)


def projected_folded(counts, c, n):
    """A variant's contribution to the folded spectrum of `n` called alleles.

    The alleles other than the major one are one synthetic minor allele, so the
    variant is treated as biallelic. A draw of `n` of the `c` called alleles
    holds `j` minor alleles with the hypergeometric chance, and `j` and `n - j`
    fall in the same bin.
    """
    major = int(np.argmax(counts))  # a tie goes to the lowest numbered allele
    minor = c - int(counts[major])
    bins = np.zeros(n // 2 + 1)
    for j in range(0, min(minor, n) + 1):
        if n - j > c - minor:
            continue
        bins[min(j, n - j)] += comb(minor, j) * comb(c - minor, n - j) / comb(c, n)
    return bins


def calc(gts, pops, ploidy, min_num_individuals, num_called_alleles):
    num_vars = gts.shape[0]
    num_alleles_seen = int(gts.max()) + 1
    g = num_called_alleles
    out = {
        name: dict(
            num_vars=0,
            num_vars_rare=0,
            num_private_vars=0,
            num_private_vars_rare=0,
            alleles=0,
            alleles_rare=0.0,
            private=0,
            private_rare=0.0,
            poly=0,
            poly_rare=0.0,
            ho=[],
            he=[],
            sfs=np.zeros((g // 2 + 1) if g else 1),
        )
        for name in pops
    }
    for v in range(num_vars):
        per_pop = {}
        for name, idx in pops.items():
            per_pop[name] = counts_of(gts[v, idx, :], num_alleles_seen)
        # A variant counts for a population when it called enough genotypes,
        # measured as the called alleles over the ploidy, as `stats` does.
        has_data = {
            name: c > 0 and c / ploidy >= min_num_individuals
            for name, (counts, c) in per_pop.items()
        }
        # The private alleles need every population to have data: an allele
        # missing from a population for want of genotypes would look private.
        all_have_data = all(has_data.values())
        all_have_draw = g is not None and all(c >= g for _, c in per_pop.values())
        for name, (counts, c) in per_pop.items():
            if not has_data[name]:
                continue
            o = out[name]
            o["num_vars"] += 1
            o["alleles"] += int((counts > 0).sum())
            if (counts > 0).sum() > 1:
                o["poly"] += 1
            gp = gts[v, pops[name], :]
            called = (gp >= 0).all(axis=1)
            if called.sum() > 0 and c > 1:
                o["ho"].append(
                    (~(gp[called] == gp[called][:, :1]).all(axis=1)).sum() / called.sum()
                )
                o["he"].append(
                    1.0 - sum(int(n) * (int(n) - 1) for n in counts) / (c * (c - 1))
                )
            if all_have_data:
                o["num_private_vars"] += 1
                others = [per_pop[q][0] for q in pops if q != name]
                o["private"] += sum(
                    1
                    for a in range(num_alleles_seen)
                    if counts[a] > 0 and all(oc[a] == 0 for oc in others)
                )
            if g is not None and c >= g:
                o["num_vars_rare"] += 1
                o["alleles_rare"] += expected_alleles(counts, c, g)
                o["poly_rare"] += chance_variable(counts, c, g)
                o["sfs"] += projected_folded(counts, c, g)
                if all_have_draw:
                    o["num_private_vars_rare"] += 1
                    others = [(per_pop[q][0], per_pop[q][1]) for q in pops if q != name]
                    o["private_rare"] += sum(
                        chance_in_draw(counts[a], c, g)
                        * np.prod(
                            [1.0 - chance_in_draw(oc[a], oc_c, g) for oc, oc_c in others]
                        )
                        for a in range(num_alleles_seen)
                        if counts[a] > 0
                    )
    res = {}
    for name, o in out.items():
        ho, he = np.array(o["ho"]), np.array(o["he"])
        nv, nr, npv = o["num_vars"], o["num_vars_rare"], o["num_private_vars"]
        npr = o["num_private_vars_rare"]
        res[name] = dict(
            num_vars=nv,
            num_vars_rare=nr,
            num_private_vars=npv,
            num_private_vars_rare=npr,
            alleles_total=o["alleles"],
            alleles_mean=o["alleles"] / nv if nv else np.nan,
            alleles_rare=o["alleles_rare"] / nr if nr else np.nan,
            private_total=o["private"],
            private_mean=o["private"] / npv if npv else np.nan,
            private_rare=o["private_rare"] / npr if npr else np.nan,
            poly_total=o["poly"],
            poly_ratio=o["poly"] / nv if nv else np.nan,
            poly_rare=o["poly_rare"] / nr if nr else np.nan,
            fis=1.0 - ho.mean() / he.mean() if len(ho) else np.nan,
            sfs=o["sfs"],
        )
    return res
