/**
 * The per variant distributions from TypeScript: the five statistics of one
 * pass, per population, and what the call refuses.
 *
 * `docs/specs/stats.md` has them under "The per variant distributions", and
 * the file they are run on is the panel of `tests/reference/stats/`, 1200
 * biallelic diploid variants of 200 individuals named `s000` to `s199`, 3 in
 * 100 genotypes missing whole, in the three populations `p0`, `p1` and `p2`
 * of 48, 68 and 84 individuals that `panel_pops_bcftools.txt` beside it
 * holds. The numbers here are the literals of `p0`, the population of
 * `s000`, that each statistic of the spec gives, and they come from the
 * `--hardy` and `--freq` reports of plink2 v2.0.0-a.7.7 on that file. The
 * comparison with pyNei itself is the one of `tests/test_stats.py`, which
 * runs both libraries; node runs neither.
 *
 * A statistic of one variant is read through the mean of a pass over that
 * one variant: `calcPerVarDistribs` gives the mean and the histogram of a
 * dataset and not the value of a variant, and the mean over one variant is
 * its value. That pass reads a VCF of the header of the panel and its first
 * data line, `var0000`, which `oneVariantOfThePanel` writes.
 */

import assert from "node:assert/strict";
import { test } from "node:test";
import { gunzipSync } from "node:zlib";

import type { PerVarDistribs, StatsDistrib, Variants } from "popnei";
import { calcPerVarDistribs, init, openVcf } from "popnei";

import { referenceStats } from "./reference.ts";

await init();

/** The bytes of the panel, gzipped, read once for the whole file. */
const PANEL = await referenceStats("panel.vcf.gz");

/** The 1200 variants of the panel. */
const PANEL_NUM_VARS = 1200;

/**
 * The individuals of each population of the panel, in the order the file of
 * populations has them: `p0` of 48 individuals, `p1` of 68 and `p2` of 84.
 */
const PANEL_POPS = popsOf(
  new TextDecoder().decode(await referenceStats("panel_pops_bcftools.txt")),
);

/** The 48 individuals of `p0`, the population of `s000`. */
const P0 = individualsOf("p0");

/** The 68 individuals of `p1`, whose unbiased value of `var0000` is known. */
const P1 = individualsOf("p1");

/**
 * The literals of `var0000` in `p0`, from the reports of plink2: 48 called
 * genotypes, 8 of them heterozygous, which is an observed heterozygosity of
 * 0.166667; 96 called alleles, 86 of the reference allele and 10 of the
 * other, which is a major allele frequency of 0.895833; and an expected
 * heterozygosity, 1 - p² - q², of 0.186632.
 */
const CALLED_GENOTYPES_OF_P0 = 48;
const OBS_HET_OF_P0 = 0.166667;
const MAF_OF_P0 = 0.895833;
const EXP_HET_OF_P0 = 0.186632;

/**
 * The unbiased expected heterozygosity of `var0000` in `p1`, which no
 * program outside the project prints: the plain 0.498998 of plink2 times
 * c / (c - 1), with c = 134 the called alleles of that population, is
 * 0.502750.
 */
const UNBIASED_EXP_HET_OF_P1 = 0.50275;

/**
 * How far a value may be from the six digits plink2 prints, one unit of the
 * last of them. The largest difference between pyNei and plink2 over the
 * 1200 variants and the 3 populations is half of it.
 */
const OF_A_PRINTED_VALUE = 1e-6;

/**
 * How far a ratio of the polymorphism counts may be from the twelve digits
 * the spec prints.
 */
const OF_A_RATIO = 1e-12;

/**
 * The polymorphism counts of the whole panel in `p0`, from the `--freq`
 * report of plink2: of its 1200 variants, all 1200 have a major allele
 * frequency in `p0`, 1173 have one below 1, which makes them variable, and
 * 1112 below 0.95, which makes them polymorphic.
 */
const POLY_OF_P0 = 1112;
const VARIABLE_OF_P0 = 1173;
const WITH_DATA_OF_P0 = 1200;
const POLY_RATIO_OF_P0 = 0.926666666667;
const POLY_RATIO_OVER_VARIABLES_OF_P0 = 0.94799658994;

/** The default histogram: 40 bins of equal width from 0 to 1, 41 edges. */
const DEFAULT_NUM_BINS = 40;

/** The individuals of each population of a file of two columns, the name of
 * an individual and the name of its population, in the order of the file. */
function popsOf(text: string): Map<string, string[]> {
  const pops = new Map<string, string[]>();
  for (const line of text.split("\n")) {
    if (line === "") {
      continue;
    }
    const [individual, pop] = line.split("\t");
    if (individual === undefined || pop === undefined) {
      throw new Error(`the line \`${line}\` is not an individual and a population`);
    }
    const individuals = pops.get(pop) ?? [];
    individuals.push(individual);
    pops.set(pop, individuals);
  }
  return pops;
}

/** The individuals of the population `pop` of the panel. */
function individualsOf(pop: string): string[] {
  const individuals = PANEL_POPS.get(pop);
  if (individuals === undefined) {
    throw new Error(`the populations of the panel hold no \`${pop}\``);
  }
  return individuals;
}

/** The 1200 variants of the panel. */
function panel(): Variants {
  return openVcf(PANEL);
}

/**
 * A VCF of `var0000` of the panel alone, the first of its data lines, with
 * the header of the panel and its 200 individuals.
 *
 * A pass over it gives one variant, so the mean of every statistic is the
 * value of that variant, which is what the literals of the spec are.
 */
function oneVariantOfThePanel(): Uint8Array {
  const lines = gunzipSync(PANEL).toString("utf8").split("\n");
  const header = lines.filter((line) => line.startsWith("#"));
  const firstVariant = lines.find((line) => line !== "" && !line.startsWith("#"));
  if (firstVariant === undefined || firstVariant.split("\t")[2] !== "var0000") {
    throw new Error("the first data line of the panel is not `var0000`");
  }
  return new TextEncoder().encode([...header, firstVariant, ""].join("\n"));
}

/** The bytes of that one variant, written once for every test. */
const VAR0000 = oneVariantOfThePanel();

/** `var0000` of the panel, over its 200 individuals. */
function theFirstVariant(): Variants {
  return openVcf(VAR0000);
}

/** The distribution of a statistic that was asked for, or the `Error` of a
 * result that holds none. */
function distribOf(
  distrib: StatsDistrib | null,
  statistic: string,
): StatsDistrib {
  if (distrib === null) {
    throw new Error(`the pass gives no ${statistic}`);
  }
  return distrib;
}

/** The value of `pop` in one array of a result, over the population names
 * `pops`. */
function ofThePop(
  values: Float64Array | Uint32Array,
  pops: readonly string[],
  pop: string,
): number {
  const value = values[pops.indexOf(pop)];
  if (value === undefined) {
    throw new Error(`the result holds no value for the population \`${pop}\``);
  }
  return value;
}

/** The mean of `statistic` in `pop`. */
function meanOf(
  distribs: PerVarDistribs,
  distrib: StatsDistrib | null,
  statistic: string,
  pop: string,
): number {
  return ofThePop(distribOf(distrib, statistic).mean, distribs.pops, pop);
}

/** Which bins of a histogram hold a count, as the number of the bin and the
 * count in it. The counts of the populations come one after another, so the
 * number of a bin of the second population is 40 and its own number. */
function binsWithACount(counts: Uint32Array): [number, number][] {
  const withACount: [number, number][] = [];
  for (const [bin, count] of counts.entries()) {
    if (count !== 0) {
      withACount.push([bin, count]);
    }
  }
  return withACount;
}

/** How many variants a histogram counted, over every bin and population. */
function variantsInTheHistogram(counts: Uint32Array): number {
  let total = 0;
  for (const count of counts) {
    total += count;
  }
  return total;
}

/** That `found` is `expected` within `tolerance`. */
function assertValue(
  found: number,
  expected: number,
  tolerance: number,
  what: string,
): void {
  assert.ok(
    Math.abs(found - expected) <= tolerance,
    `${what} is ${found} and not ${expected} within ${tolerance}`,
  );
}

test("the statistics of var0000 in p0 are the literals of plink2", () => {
  // The pass gives one variant, so the mean of each statistic is the value
  // of that variant in that population. The unbiased expected
  // heterozygosity is read in `p1`, the one population a number is known
  // for: no program outside the project prints it.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, { pops: { p0: P0, p1: P1 } });

  assert.deepEqual(distribs.pops, ["p0", "p1"]);
  assertValue(
    meanOf(distribs, distribs.obsHet, "observed heterozygosity", "p0"),
    OBS_HET_OF_P0,
    OF_A_PRINTED_VALUE,
    "the observed heterozygosity of var0000 in p0",
  );
  assertValue(
    meanOf(distribs, distribs.maf, "major allele frequency", "p0"),
    MAF_OF_P0,
    OF_A_PRINTED_VALUE,
    "the maf of var0000 in p0",
  );
  assertValue(
    meanOf(distribs, distribs.expHet, "expected heterozygosity", "p0"),
    EXP_HET_OF_P0,
    OF_A_PRINTED_VALUE,
    "the plain expected heterozygosity of var0000 in p0",
  );
  assertValue(
    meanOf(
      distribs,
      distribs.unbiasedExpHet,
      "unbiased expected heterozygosity",
      "p1",
    ),
    UNBIASED_EXP_HET_OF_P1,
    OF_A_PRINTED_VALUE,
    "the unbiased expected heterozygosity of var0000 in p1",
  );
  assert.deepEqual(distribs.passStats, { numVars: 1, filtering: {} });
  variants.free();
});

test("the histogram of the pass counts each value in its bin, one population after another", () => {
  // The 40 bins of the default histogram are 0.025 wide, so the observed
  // heterozygosity of 0.166667 of `p0` falls in the seventh of them, 0.15
  // to 0.175, which is the bin 6, and the 0.507463 of `p1` in the
  // twenty-first, 0.5 to 0.525, which is the bin 20. The counts of `p1`
  // come after the 40 of `p0`, so its bin is the 60th number of the array.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, {
    stats: ["obs_het"],
    pops: { p0: P0, p1: P1 },
  });

  const obsHet = distribOf(distribs.obsHet, "observed heterozygosity");
  assert.equal(obsHet.histBinEdges.length, DEFAULT_NUM_BINS + 1);
  assert.equal(obsHet.histBinEdges[0], 0);
  assert.equal(obsHet.histBinEdges[DEFAULT_NUM_BINS], 1);
  assert.equal(obsHet.histCounts.length, 2 * DEFAULT_NUM_BINS);
  assert.deepEqual(binsWithACount(obsHet.histCounts), [
    [6, 1],
    [60, 1],
  ]);
  variants.free();
});

test("the polymorphism counts of the panel in p0 are plink2's", () => {
  // Of the 1200 variants, 1112 have a major allele frequency below 0.95 in
  // `p0` and 1173 one below 1, and all 1200 have one: every population of
  // the panel has 40 called genotypes or more at every variant, above the
  // 20 of the default threshold, so every variant has a value in each of
  // them and every histogram counts the 1200.
  const variants = panel();

  const distribs = calcPerVarDistribs(variants, { pops: { p0: P0 } });

  const poly = distribs.polyVarsRatio;
  if (poly === null) {
    throw new Error("the pass gives no polymorphism counts");
  }
  assert.deepEqual(poly.numPoly, Uint32Array.of(POLY_OF_P0));
  assert.deepEqual(poly.numVariable, Uint32Array.of(VARIABLE_OF_P0));
  assert.deepEqual(poly.totNumVariantsWithData, Uint32Array.of(WITH_DATA_OF_P0));
  assertValue(
    ofThePop(poly.polyRatio, distribs.pops, "p0"),
    POLY_RATIO_OF_P0,
    OF_A_RATIO,
    "the polymorphism ratio of the panel in p0",
  );
  assertValue(
    ofThePop(poly.polyRatioOverVariables, distribs.pops, "p0"),
    POLY_RATIO_OVER_VARIABLES_OF_P0,
    OF_A_RATIO,
    "the polymorphism ratio over the variable variants of the panel in p0",
  );
  assert.equal(distribs.passStats.numVars, PANEL_NUM_VARS);
  for (const [distrib, statistic] of [
    [distribs.obsHet, "observed heterozygosity"],
    [distribs.maf, "major allele frequency"],
    [distribs.expHet, "expected heterozygosity"],
    [distribs.unbiasedExpHet, "unbiased expected heterozygosity"],
  ] as const) {
    assert.equal(
      variantsInTheHistogram(distribOf(distrib, statistic).histCounts),
      PANEL_NUM_VARS,
      `the histogram of the ${statistic} counts every variant of the panel`,
    );
  }
  variants.free();
});

test("a pass with no pops has one population of the individuals its steps give", () => {
  // The names of the populations are resolved against the individuals the
  // pass gives, which are those a filter of individuals kept, and with no
  // `pops` the one population is those individuals: the counts of the 48
  // individuals of `p0` kept by the filter are the counts of `p0`.
  const variants = panel();
  variants.filterIndividuals(P0);

  const distribs = calcPerVarDistribs(variants, {
    stats: ["poly_vars_ratio"],
  });

  assert.deepEqual(distribs.pops, ["pop"]);
  const poly = distribs.polyVarsRatio;
  if (poly === null) {
    throw new Error("the pass gives no polymorphism counts");
  }
  assert.deepEqual(poly.numPoly, Uint32Array.of(POLY_OF_P0));
  assert.deepEqual(poly.numVariable, Uint32Array.of(VARIABLE_OF_P0));
  assert.deepEqual(poly.totNumVariantsWithData, Uint32Array.of(WITH_DATA_OF_P0));
  variants.free();
});

test("a population of 15 individuals has no value at the default threshold", () => {
  // `minNumIndividuals` is 20 by default, so a population of 15 has fewer
  // called genotypes than that at every variant: every mean is NaN, every
  // histogram is empty and the polymorphism counts are 0, with both ratios
  // NaN.
  const variants = panel();

  const distribs = calcPerVarDistribs(variants, {
    pops: { small: P0.slice(0, 15) },
  });

  for (const [distrib, statistic] of [
    [distribs.obsHet, "observed heterozygosity"],
    [distribs.maf, "major allele frequency"],
    [distribs.expHet, "expected heterozygosity"],
    [distribs.unbiasedExpHet, "unbiased expected heterozygosity"],
  ] as const) {
    const ofTheStat = distribOf(distrib, statistic);
    assert.ok(
      Number.isNaN(ofTheStat.mean[0]),
      `the mean of the ${statistic} of a population of 15 is not NaN`,
    );
    assert.equal(variantsInTheHistogram(ofTheStat.histCounts), 0);
  }
  const poly = distribs.polyVarsRatio;
  if (poly === null) {
    throw new Error("the pass gives no polymorphism counts");
  }
  assert.deepEqual(poly.numPoly, Uint32Array.of(0));
  assert.deepEqual(poly.numVariable, Uint32Array.of(0));
  assert.deepEqual(poly.totNumVariantsWithData, Uint32Array.of(0));
  assert.ok(Number.isNaN(poly.polyRatio[0]));
  assert.ok(Number.isNaN(poly.polyRatioOverVariables[0]));
  variants.free();
});

test("a variant with fewer called genotypes than minNumIndividuals has no value", () => {
  // `p0` has 48 called genotypes at `var0000`, and a variant has no value
  // when its called data, counted in genotypes, is strictly below the
  // threshold: at 48 the statistics are those of plink2 and at 49 there are
  // none.
  const kept = theFirstVariant();
  const dropped = theFirstVariant();

  const atTheThreshold = calcPerVarDistribs(kept, {
    pops: { p0: P0 },
    minNumIndividuals: CALLED_GENOTYPES_OF_P0,
  });
  const aboveIt = calcPerVarDistribs(dropped, {
    pops: { p0: P0 },
    minNumIndividuals: CALLED_GENOTYPES_OF_P0 + 1,
  });

  assertValue(
    meanOf(atTheThreshold, atTheThreshold.obsHet, "observed heterozygosity", "p0"),
    OBS_HET_OF_P0,
    OF_A_PRINTED_VALUE,
    "the observed heterozygosity of var0000 in p0 at a threshold of 48",
  );
  const withoutAValue = distribOf(aboveIt.obsHet, "observed heterozygosity");
  assert.ok(Number.isNaN(withoutAValue.mean[0]));
  assert.equal(variantsInTheHistogram(withoutAValue.histCounts), 0);
  kept.free();
  dropped.free();
});

test("only the statistics that were asked for are calculated", () => {
  // Asking for fewer is a saving of work and changes no value, and what
  // nobody asked for is null.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, {
    stats: ["maf", "poly_vars_ratio"],
    pops: { p0: P0 },
  });

  assert.equal(distribs.obsHet, null);
  assert.equal(distribs.expHet, null);
  assert.equal(distribs.unbiasedExpHet, null);
  assert.notEqual(distribs.polyVarsRatio, null);
  assertValue(
    meanOf(distribs, distribs.maf, "major allele frequency", "p0"),
    MAF_OF_P0,
    OF_A_PRINTED_VALUE,
    "the maf of var0000 in p0",
  );
  variants.free();
});

test("an exponent of 1 leaves both expected heterozygosities at 0", () => {
  // `ploidy` is the exponent of the frequencies and the number of factors
  // of the products of the unbiased one. At 1 the plain one is 1 less the
  // sum of the frequencies and the unbiased one is 1 less the sum of the
  // counts over the called alleles, and both sums are 1 at every variant
  // with a called allele.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, {
    pops: { p0: P0 },
    ploidy: 1,
  });

  assertValue(
    meanOf(distribs, distribs.expHet, "expected heterozygosity", "p0"),
    0,
    OF_A_RATIO,
    "the plain expected heterozygosity of var0000 in p0 at an exponent of 1",
  );
  assertValue(
    meanOf(
      distribs,
      distribs.unbiasedExpHet,
      "unbiased expected heterozygosity",
      "p0",
    ),
    0,
    OF_A_RATIO,
    "the unbiased expected heterozygosity of var0000 in p0 at an exponent of 1",
  );
  variants.free();
});

test("logarithmic bins of equal ratio span the range they were given", () => {
  // The example of `test_hist.py` of pyNei: a logarithmic range from 0.01
  // to 100 in 4 bins has the edges 0.01, 0.1, 1, 10 and 100. Two libraries
  // need not round the powers of 10 alike, so the edges are compared within
  // 1e-12.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, {
    stats: ["maf"],
    pops: { p0: P0 },
    histKwargs: { range: [0.01, 100], numBins: 4, binType: "logarithmic" },
  });

  const edges = distribOf(distribs.maf, "major allele frequency").histBinEdges;
  assert.equal(edges.length, 5);
  for (const [edge, expected] of [...edges].map(
    (found, bin) => [found, [0.01, 0.1, 1, 10, 100][bin]] as const,
  )) {
    assertValue(edge, expected ?? Number.NaN, OF_A_RATIO, "an edge of the bins");
  }
  variants.free();
});

test("a poly threshold below the frequency of a variant leaves it out of the polymorphic ones", () => {
  // `var0000` has a major allele frequency of 0.895833 in `p0`, so it is
  // polymorphic below the default threshold of 0.95 and not below 0.5. It
  // is variable either way, since that frequency is below 1.
  const variants = theFirstVariant();

  const distribs = calcPerVarDistribs(variants, {
    stats: ["poly_vars_ratio"],
    pops: { p0: P0 },
    polyThreshold: 0.5,
  });

  const poly = distribs.polyVarsRatio;
  if (poly === null) {
    throw new Error("the pass gives no polymorphism counts");
  }
  assert.deepEqual(poly.numPoly, Uint32Array.of(0));
  assert.deepEqual(poly.numVariable, Uint32Array.of(1));
  assert.deepEqual(poly.totNumVariantsWithData, Uint32Array.of(1));
  assert.equal(poly.polyRatio[0], 0);
  assert.equal(poly.polyRatioOverVariables[0], 0);
  variants.free();
});

test("a name that is not an individual of the pass is refused", () => {
  // The names are looked up among the individuals the pass gives, and one
  // that is none of them is an `Error` that names it: a population of
  // whoever was found would be a result nobody asked for.
  const variants = panel();

  assert.throws(
    () => calcPerVarDistribs(variants, { pops: { p0: ["s000", "nobody"] } }),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("nobody") &&
      error.message.includes("p0"),
  );

  variants.free();
});

test("a key of histKwargs that popnei does not know is refused", () => {
  // pyNei ignores it, so a user who writes `nbins` gets the 40 bins of the
  // default with nothing said. The message names the three keys there are.
  const variants = theFirstVariant();

  assert.throws(
    () =>
      calcPerVarDistribs(variants, {
        pops: { p0: P0 },
        histKwargs: { nbins: 4 } as unknown as { numBins?: number },
      }),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("nbins") &&
      error.message.includes("numBins"),
  );

  variants.free();
});

test("a kind of bins that is neither of the two is refused under its name", () => {
  // The core knows the two kinds and names the argument `bin_type`, which
  // is what a Python user writes; a TypeScript user wrote `binType`, and
  // that is the call they have to look at. `lineal` is how pyNei spells the
  // bins of equal width, and popnei refuses it as any other unknown name.
  const variants = theFirstVariant();

  assert.throws(
    () =>
      calcPerVarDistribs(variants, {
        pops: { p0: P0 },
        histKwargs: { binType: "lineal" as unknown as "linear" },
      }),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("`binType` is `lineal`") &&
      error.message.includes("logarithmic") &&
      !error.message.includes("bin_type"),
  );

  variants.free();
});

test("a pass that calculates no statistic at all is refused", () => {
  // A result holds the statistics that were asked for, and one that holds
  // none is a pass over the whole file for nothing.
  const variants = theFirstVariant();

  assert.throws(
    () => calcPerVarDistribs(variants, { stats: [], pops: { p0: P0 } }),
    (error: unknown) => error instanceof Error && error.message.includes("stats"),
  );

  variants.free();
});

test("a statistic that is not one of the five is refused", () => {
  // The five are a union of string literals in TypeScript, so a typo does
  // not compile; what reaches the call is a name written in JavaScript.
  const variants = theFirstVariant();

  assert.throws(
    () =>
      calcPerVarDistribs(variants, {
        stats: ["obs_hets"] as unknown as ["obs_het"],
        pops: { p0: P0 },
      }),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("obs_hets") &&
      error.message.includes("poly_vars_ratio"),
  );

  variants.free();
});
