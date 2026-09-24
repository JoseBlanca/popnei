/**
 * How much variety each population holds, from TypeScript: the alleles each
 * population called, the private ones among them, the variants that vary in
 * it, how its variants are spread over the count of their rarer allele and
 * F_IS, over the panel and over the cases a user can reach.
 *
 * `docs/specs/diversity.md` has the five statistics and, under "How it is
 * verified" of each, the program its numbers come from. The literals here
 * are the ones `tests/test_diversity.py` asserts, which is what goal 1 of
 * `docs/objectives.md` asks for: the two packages give a user the same
 * numbers. The Python suite reads them from the files of
 * `tests/reference/diversity/`, which `adegenet` 2.1.11, `poppr` 2.9.8,
 * `vegan` 2.7.6, `dadi` 2.4.4 and `scikit-allel` 1.3.13 wrote; node runs none
 * of those programs, so they are written here as literals, as section 11 of
 * `docs/architecture.md` has it for the tests of this package.
 *
 * The dataset is the panel of `docs/specs/stats.md`,
 * `tests/reference/stats/panel.vcf.gz`: 1200 biallelic diploid variants of
 * 200 individuals, 3 in 100 genotypes missing whole, in the three
 * populations `p0`, `p1` and `p2` of 48, 68 and 84 individuals that
 * `panel_pops_bcftools.txt` beside it holds. At `minNumIndividuals` 20 all
 * 1200 variants count for all three.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PopDiversity, Variants } from "popnei";
import {
  calcPopDiversity,
  init,
  openVcf,
  popDiversityStatsWithoutADraw,
} from "popnei";

import { referenceStats, vcfOf } from "./reference.ts";

await init();

/** The bytes of the panel, gzipped, read once for the whole file. */
const PANEL = await referenceStats("panel.vcf.gz");

/** The 1200 variants of the panel. */
const PANEL_NUM_VARS = 1200;

/**
 * How many called genotypes a population needs at a variant for the variant
 * to count for it, which is what "How it is verified" of every item of the
 * spec ran the reference programs at.
 */
const PANEL_MIN_NUM_INDIVIDUALS = 20;

/**
 * The individuals of each population of the panel, under the name of the
 * population, in the order `panel_pops_bcftools.txt` names them, which is the
 * order every result of a call with these populations is in.
 */
const PANEL_POPS = popsOf(
  new TextDecoder().decode(await referenceStats("panel_pops_bcftools.txt")),
);

/**
 * The three populations in the order of the keys of `PANEL_POPS`, which is
 * the order `panel_pops_bcftools.txt` names them in and the order every
 * result of a call with those populations has to be in.
 *
 * It is not alphabetical, and that is what it is for: `p2` comes before `p1`
 * in the file, so every array of literals below is written in that order, and
 * a result that put the populations in any other would fail here. The Python
 * suite reads the same file and keeps the same order.
 */
const PANEL_POP_NAMES = ["p0", "p2", "p1"];

/**
 * The alleles the three populations called over the 1200 variants and the
 * mean of each, which is the allelic richness, from `adegenet`.
 */
const PANEL_NUM_ALLELES = [2373, 2384, 2377];
const PANEL_MEAN_NUM_ALLELES = [
  1.9775, 1.9866666666666666, 1.9808333333333332,
];

/**
 * Of those, the ones no other population of the call called at the same
 * variant, from `poppr`, and their mean over the 1200 variants that counted
 * for every population.
 */
const PANEL_PRIVATE_ALLELES = [0, 1, 0];
const PANEL_MEAN_PRIVATE_ALLELES = [0, 0.0008333333333333334, 0];

/**
 * The variants each population called more than one allele at, from
 * `adegenet`, and that count over the 1200 variants that counted for the
 * population.
 */
const PANEL_VARIABLE_VARS = [1173, 1184, 1177];
const PANEL_VARIABLE_VARS_RATIO = [
  0.9775, 0.9866666666666667, 0.9808333333333333,
];

/**
 * The unbiased F_IS of the three populations, from "How it is verified" of
 * "The inbreeding coefficient F_IS" of the spec. It is the form
 * `calcPopDiversity` gives: one minus the mean observed heterozygosity over
 * the mean unbiased expected one.
 *
 * The item prints those three to ten decimals for a reader and gives them
 * again at the precision they were computed to, which is what is written
 * here: a value of ten decimals stands for anything within 5e-11 of itself,
 * and 1e-12 of 0.0128 is 1.3e-14, so the shorter form could not be compared
 * within the tolerance this file uses. They come from
 * `docs/reports/diversity-method/panel.py`, which computes the five
 * quantities of the spec in Python, so this assertion says that popnei
 * agrees with that Python and nothing more. What checks F_IS against a
 * program outside popnei is the plain form of it, which
 * `tests/test_diversity.py` builds from `calcPerVarDistribs` and compares
 * with `scikit-allel`.
 */
const PANEL_FIS = [
  -0.012758486763377208, -0.018458583231322434, -0.018110713076467277,
];

/**
 * How far a value of popnei may be from the number it is compared with,
 * which is what every item of the spec asks for its floats: 1e-12 of the
 * value. Every number here carries all its digits, and the two sides add the
 * same per variant values in different orders, so the last bits differ.
 */
const OF_A_REFERENCE_VALUE = 1e-12;

/**
 * How many called alleles every population is brought down to in the draw the
 * reference programs were run at, and how many bins a spectrum of that draw
 * has: the counts of the rarer allele from 0 to 20 / 2.
 */
const PANEL_NUM_CALLED_ALLELES = 20;
const PANEL_SFS_BINS = 11;

/**
 * Every gene copy the panel holds, its 200 individuals at a ploidy of 2,
 * which is the largest draw it allows: the 3 in 100 genotypes it is missing
 * leave no population able to call that many alleles at any variant, so the
 * draw is taken and every value of it is missing, and one allele more is
 * refused.
 */
const EVERY_GENE_COPY_OF_THE_PANEL = 400;

/**
 * The alleles a draw of 20 is expected to show in each population, which
 * `vegan` measured, and the chance that such a draw shows more than one
 * allele, which is each of those less 1: every variant of the panel has two
 * alleles, so a draw there shows one of them or both.
 */
const PANEL_NUM_ALLELES_IN_DRAW = [
  1.9283948650041205, 1.9197370843937562, 1.9219209943237829,
];
const PANEL_VARIABLE_VARS_RATIO_IN_DRAW = PANEL_NUM_ALLELES_IN_DRAW.map(
  (alleles) => alleles - 1,
);

/**
 * The alleles a draw of 20 is expected to show in one population and in no
 * other, averaged over the 1200 variants every population reached the draw
 * at, from "How it is verified" of "The private alleles" of the spec, which
 * gives them to ten decimals.
 *
 * No program outside popnei computes a standardized private allele value, so
 * these come from `docs/reports/diversity-method/panel.py` as the unbiased
 * F_IS above does.
 */
const PANEL_PRIVATE_ALLELES_IN_DRAW = [
  0.0112196177, 0.0089014974, 0.0099715392,
];

/**
 * What a number printed to ten decimals stands for, which is the bound those
 * three are compared within: anything within 5e-11 of what is written rounds
 * to the same ten decimals.
 */
const OF_TEN_DECIMALS = 5e-11;

/**
 * How many variants of each population a draw of 20 is expected to show each
 * count of the rarer allele at, from `dadi` through
 * `tests/reference/diversity/panel_folded_sfs_dadi.tsv`: one array per
 * population, the count 0 first.
 *
 * Each column sums to the 1200 variants that counted, and the sum is compared
 * within the tolerance and not exactly: `dadi`'s stored columns are short of
 * 1200 by 1.3e-11, 7.0e-11 and 4.0e-11 and popnei's own are within 2.3e-13,
 * both inside the 1.2e-9 that 1e-12 of 1200 allows, measured on 24 September
 * 2026.
 */
const PANEL_FOLDED_SFS = [
  [
    85.92616199505309, 92.99651940325519, 106.88963282305795,
    115.54543767638108, 120.50559386225342, 122.94119349597617,
    124.08151270035586, 124.23937415023605, 123.49607649027779,
    122.42162180189403, 60.95687560124602,
  ],
  [
    96.3154987274892, 101.37814416196603, 108.05152927646787,
    114.8500079366327, 119.73977861754621, 121.88331138362474,
    121.86450896763563, 120.60097546763699, 118.97681357034003,
    117.71473832544012, 58.62469356518016,
  ],
  [
    93.69480681145635, 95.45880325926598, 103.72623430085372,
    110.6831031478722, 116.6729455230275, 121.05177240834107,
    123.60219049366185, 124.59443950850212, 124.54032617876919,
    124.0682314877391, 61.907146880441395,
  ],
];

/**
 * The four statistics that need no draw, which are every one but the folded
 * spectrum: the bins of a spectrum are the counts of the rarer allele in a
 * draw, so a call that names it gives `numCalledAlleles` as well.
 */
const WITH_NO_SPECTRUM = [
  "num_alleles",
  "private_alleles",
  "variable_vars_ratio",
  "fis",
] as const;

/**
 * Four variants of the three diploid individuals `ind1`, `ind2` and `ind3`,
 * for the counts that the panel cannot show: there every population counted
 * every variant, so the variants of a population and the variants of every
 * population are the same 1200 and a divisor read from the wrong one of the
 * two gives the same number.
 *
 * The genotypes are `0/0 1/1 1/1`, `0/0 0/0 0/0`, `0/1 0/0 0/0` and
 * `0/0 ./. ./.`: the last one is called in `ind1` alone, so a population of
 * `ind2` and `ind3` counts three variants where a population of `ind1`
 * counts four.
 */
const SOME_VARS_MISSED = vcfOf([
  "chr1\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t1/1",
  "chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0",
  "chr1\t3\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0",
  "chr1\t4\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t./.\t./.",
]);

/**
 * Two variants of the same three individuals, for the four counts of variants
 * that the panel cannot tell apart either: there all four are 1200, so a
 * divisor or a count read from the wrong one of them gives the same number.
 *
 * The genotypes are `0/1 1/1 2/2` and `0/0 0/1 0/.`. With `ind1` and `ind2` in
 * one population and `ind3` in another, at a draw of 2 and a
 * `minNumIndividuals` of 0, the first population calls 4 alleles at both
 * variants and the second calls 2 at the first and, its genotype there being
 * the half called `0/.`, one at the second: that one allele counts the variant
 * for it, the population having called something, and is below the draw. So
 * the variants with data are 2 and 2, the variants in the draw 2 and 1, the
 * variants every population counted 2 and the variants every population
 * reached the draw at 1.
 */
const SOME_VARS_SHORT_OF_THE_DRAW = vcfOf([
  "chr1\t10\t.\tA\tT,G\t.\tPASS\t.\tGT\t0/1\t1/1\t2/2",
  "chr1\t20\t.\tA\tT,G\t.\tPASS\t.\tGT\t0/0\t0/1\t0/.",
]);

/**
 * The individuals of each population of a file of two columns, the name of
 * an individual and the name of its population, under the names of the
 * populations in the order the file names them.
 *
 * They are not sorted. The order of the keys is the order of every array of a
 * result, and the file names `p2` before `p1`, so a result that sorted the
 * populations anywhere, or that put one population's counts under another's
 * name, gives the literals above in the wrong order and fails.
 */
function popsOf(text: string): Record<string, string[]> {
  const ofEachName = new Map<string, string[]>();
  for (const line of text.split("\n")) {
    if (line === "") {
      continue;
    }
    const [individual, pop] = line.split("\t");
    if (individual === undefined || pop === undefined) {
      throw new Error(
        `the line \`${line}\` is not an individual and a population`,
      );
    }
    const individuals = ofEachName.get(pop) ?? [];
    individuals.push(individual);
    ofEachName.set(pop, individuals);
  }
  return Object.fromEntries(ofEachName);
}

/** The 1200 variants of the panel. */
function panel(): Variants {
  return openVcf(PANEL);
}

/**
 * The diversity of the three populations of the panel, with the threshold
 * the reference programs were run at, and the variants freed when it
 * returns.
 */
function ofThePanel(
  options: Parameters<typeof calcPopDiversity>[1] = {},
): PopDiversity {
  const variants = panel();
  try {
    return calcPopDiversity(variants, {
      pops: PANEL_POPS,
      minNumIndividuals: PANEL_MIN_NUM_INDIVIDUALS,
      ...options,
    });
  } finally {
    variants.free();
  }
}

/** That `found` holds the numbers of `expected`, each within
 * `OF_A_REFERENCE_VALUE` of its value. */
function assertValues(
  found: Float64Array | Uint32Array | null | undefined,
  expected: readonly number[],
  what: string,
): void {
  if (found === null || found === undefined) {
    throw new Error(`the result holds no ${what}`);
  }
  assert.equal(found.length, expected.length, `the populations of ${what}`);
  for (const [pop, value] of expected.entries()) {
    const ours = found[pop] as number;
    assert.ok(
      Math.abs(ours - value) <=
        OF_A_REFERENCE_VALUE * Math.max(Math.abs(ours), Math.abs(value)),
      `the ${what} of ${PANEL_POP_NAMES[pop]} is ${ours} and not ${value}`,
    );
  }
}

/**
 * That `found` holds the numbers of `expected`, each within `OF_TEN_DECIMALS`
 * of its value, which is what a number the spec prints to ten decimals stands
 * for.
 */
function assertToTenDecimals(
  found: Float64Array,
  expected: readonly number[],
  what: string,
): void {
  assert.equal(found.length, expected.length, `the populations of ${what}`);
  for (const [pop, value] of expected.entries()) {
    const ours = found[pop] as number;
    assert.ok(
      Math.abs(ours - value) <= OF_TEN_DECIMALS,
      `the ${what} of ${PANEL_POP_NAMES[pop]} is ${ours} and not ${value}`,
    );
  }
}

/** The counts of a statistic that was asked for, or the `Error` of a result
 * that holds none. */
function counted<T>(counts: T | null, statistic: string): T {
  if (counts === null) {
    throw new Error(`the pass gives no ${statistic}`);
  }
  return counts;
}

test("the alleles of the panel and its F_IS are the numbers of the reference programs", () => {
  // All 1200 variants of the panel count for all three populations at a
  // threshold of 20 called genotypes, so the counts below are over the
  // whole panel and the private alleles are over the same 1200 variants.
  const diversity = ofThePanel({ stats: WITH_NO_SPECTRUM });

  assert.deepEqual(diversity.pops, PANEL_POP_NAMES);
  assert.deepEqual(
    counted(diversity.numAlleles, "count of alleles").total,
    Uint32Array.from(PANEL_NUM_ALLELES),
  );
  assert.deepEqual(
    counted(diversity.privateAlleles, "count of private alleles").total,
    Uint32Array.from(PANEL_PRIVATE_ALLELES),
  );
  assert.deepEqual(
    counted(diversity.variableVarsRatio, "count of variable variants").total,
    Uint32Array.from(PANEL_VARIABLE_VARS),
  );
  assertValues(diversity.fis, PANEL_FIS, "F_IS");
  assert.deepEqual(
    diversity.numVars.withData,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
  assert.equal(diversity.numVarsEveryPop, PANEL_NUM_VARS);
  assert.equal(diversity.passStats.numVars, PANEL_NUM_VARS);
  assert.deepEqual(diversity.passStats.filtering, {});
});

test("the means and the ratios of the panel are the totals over the variants that counted", () => {
  // The mean alleles and the ratio of variable variants are over the
  // variants that counted for the population, and the mean private alleles
  // over the variants that counted for every population.
  const diversity = ofThePanel({ stats: WITH_NO_SPECTRUM });

  assertValues(
    counted(diversity.numAlleles, "count of alleles").mean,
    PANEL_MEAN_NUM_ALLELES,
    "mean alleles",
  );
  assertValues(
    counted(diversity.privateAlleles, "count of private alleles").mean,
    PANEL_MEAN_PRIVATE_ALLELES,
    "mean private alleles",
  );
  assertValues(
    counted(diversity.variableVarsRatio, "count of variable variants").ratio,
    PANEL_VARIABLE_VARS_RATIO,
    "ratio of variable variants",
  );
});

test("a statistic that was not asked for has no value", () => {
  // The counts of the variants are there whatever was asked for: they are
  // what every mean and every ratio is over.
  const diversity = ofThePanel({ stats: ["fis"] });

  assertValues(diversity.fis, PANEL_FIS, "F_IS");
  assert.equal(diversity.numAlleles, null);
  assert.equal(diversity.privateAlleles, null);
  assert.equal(diversity.variableVarsRatio, null);
  assert.equal(diversity.foldedSfs, null);
  assert.deepEqual(
    diversity.numVars.withData,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
});

test("one population holds every allele it called as a private one", () => {
  // With no `pops` there is one population of every individual, and with no
  // other population to hold them every allele it called is private.
  const diversity = ofThePanel({
    pops: undefined,
    stats: ["num_alleles", "private_alleles"],
  });

  assert.deepEqual(diversity.pops, ["pop"]);
  const total = counted(diversity.numAlleles, "count of alleles").total;
  assert.ok((total[0] as number) > 0);
  assert.deepEqual(
    counted(diversity.privateAlleles, "count of private alleles").total,
    total,
  );
  assert.equal(diversity.numVarsEveryPop, PANEL_NUM_VARS);
});

test("a call that names no statistic gives the four that need no draw", () => {
  // `calcPopDiversity(variants)` with no options at all: the four statistics
  // that need no draw and no spectrum, over one population of every
  // individual. The bins of the spectrum are counts of the rarer allele in a
  // draw of `numCalledAlleles`, so while the default was all five this call
  // refused itself and named an option the user had not written.
  const variants = panel();
  let diversity: PopDiversity;
  try {
    diversity = calcPopDiversity(variants);
  } finally {
    variants.free();
  }

  assert.deepEqual(diversity.pops, ["pop"]);
  const alleles = counted(diversity.numAlleles, "count of alleles");
  assert.ok((alleles.total[0] as number) > 0);
  assert.ok(
    (counted(diversity.privateAlleles, "count of private alleles")
      .total[0] as number) > 0,
  );
  assert.ok(
    (counted(diversity.variableVarsRatio, "ratio of variable variants")
      .total[0] as number) > 0,
  );
  assert.ok(!Number.isNaN((diversity.fis as Float64Array)[0] as number));
  // The one statistic that needs a draw is not in the default, and naming it
  // without a draw is refused as it was.
  assert.equal(diversity.foldedSfs, null);
  assert.deepEqual(diversity.numVars.withData, Uint32Array.of(PANEL_NUM_VARS));
});

test("the four statistics that need no draw are the default and have a name", () => {
  // The four are named in the core, and a user who wants them and the folded
  // spectrum takes them from here rather than writing the names into their own
  // code, as a Python user takes `PopDiversityStat.WITHOUT_A_DRAW`.
  assert.deepEqual(popDiversityStatsWithoutADraw(), [...WITH_NO_SPECTRUM]);

  const diversity = ofThePanel({
    stats: [...popDiversityStatsWithoutADraw(), "folded_sfs"],
    numCalledAlleles: PANEL_NUM_CALLED_ALLELES,
  });

  assertValues(
    counted(diversity.numAlleles, "count of alleles").inDraw,
    PANEL_NUM_ALLELES_IN_DRAW,
    "alleles a draw of 20 shows",
  );
  assertValues(diversity.fis, PANEL_FIS, "F_IS");
  assert.equal(
    Object.keys(counted(diversity.foldedSfs, "folded spectrum")).length,
    PANEL_POP_NAMES.length,
  );
});

test("the standardized values of the panel in a draw of 20 are the ones of the reference", () => {
  // The alleles a draw of 20 shows and the chance that such a draw varies,
  // which `vegan` measured, and the private alleles of the draw, which the
  // spec gives: all three populations of the panel reach 20 called alleles at
  // every one of its 1200 variants, so the three standardized values are
  // means over all of them and `numVars.inDraw` is 1200.
  const diversity = ofThePanel({
    stats: WITH_NO_SPECTRUM,
    numCalledAlleles: PANEL_NUM_CALLED_ALLELES,
  });

  assertValues(
    counted(diversity.numAlleles, "count of alleles").inDraw,
    PANEL_NUM_ALLELES_IN_DRAW,
    "alleles a draw of 20 shows",
  );
  assertValues(
    counted(diversity.variableVarsRatio, "count of variable variants").inDraw,
    PANEL_VARIABLE_VARS_RATIO_IN_DRAW,
    "chance that a draw of 20 varies",
  );
  assertToTenDecimals(
    counted(diversity.privateAlleles, "count of private alleles").inDraw,
    PANEL_PRIVATE_ALLELES_IN_DRAW,
    "private alleles a draw of 20 shows",
  );
  assert.deepEqual(
    diversity.numVars.inDraw,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
  assert.equal(diversity.numVarsEveryPopInDraw, PANEL_NUM_VARS);
  // A draw changes the standardized values alone: the totals and F_IS are the
  // numbers of a call that gave none.
  const withNoDraw = ofThePanel({ stats: WITH_NO_SPECTRUM });
  assert.deepEqual(
    counted(diversity.numAlleles, "count of alleles").total,
    counted(withNoDraw.numAlleles, "count of alleles").total,
  );
  assert.deepEqual(
    counted(diversity.numAlleles, "count of alleles").mean,
    counted(withNoDraw.numAlleles, "count of alleles").mean,
  );
  assert.deepEqual(diversity.fis, withNoDraw.fis);
});

test("the folded spectrum of the panel in a draw of 20 is the one dadi projected", () => {
  // Eleven bins, the counts of the rarer allele 0 to 10, one `Float64Array`
  // per population under its name. Each variant in the draw gives every bin
  // the chance that a draw of 20 shows that many rarer copies there, so the
  // values are not whole numbers and each column sums to the 1200 variants it
  // was taken over.
  const diversity = ofThePanel({
    stats: ["folded_sfs"],
    numCalledAlleles: PANEL_NUM_CALLED_ALLELES,
  });

  const spectra = counted(diversity.foldedSfs, "folded spectrum");
  assert.deepEqual(Object.keys(spectra), PANEL_POP_NAMES);
  for (const [pop, name] of PANEL_POP_NAMES.entries()) {
    const ours = spectra[name] as Float64Array;
    const theirs = PANEL_FOLDED_SFS[pop] as number[];
    assert.equal(ours.length, PANEL_SFS_BINS, `the bins of ${name}`);
    for (const [rarerAllele, value] of theirs.entries()) {
      const found = ours[rarerAllele] as number;
      assert.ok(
        Math.abs(found - value) <=
          OF_A_REFERENCE_VALUE * Math.max(Math.abs(found), Math.abs(value)),
        `the variants of ${name} with ${rarerAllele} copies of the rarer ` +
          `allele are ${found} and not ${value}`,
      );
    }
    const total = ours.reduce((sum, value) => sum + value, 0);
    assert.ok(
      Math.abs(total - PANEL_NUM_VARS) <=
        OF_A_REFERENCE_VALUE * PANEL_NUM_VARS,
      `the spectrum of ${name} sums to ${total} and not ${PANEL_NUM_VARS}`,
    );
  }
  // The other four statistics were not asked for, and the counts of the
  // variants are there whatever was asked for.
  assert.equal(diversity.numAlleles, null);
  assert.deepEqual(
    diversity.numVars.inDraw,
    Uint32Array.of(PANEL_NUM_VARS, PANEL_NUM_VARS, PANEL_NUM_VARS),
  );
});

test("the spectrum asked for without a draw is refused", () => {
  // The bins of a folded spectrum are the counts of the rarer allele in a
  // draw, so the spectrum needs `numCalledAlleles`. The message is the core's,
  // with the name of the option written as a TypeScript user wrote it.
  assert.throws(() => ofThePanel({ stats: ["folded_sfs"] }), {
    message: /`numCalledAlleles` was not given/,
  });
  assert.throws(() => ofThePanel({ stats: ["fis", "folded_sfs"] }), {
    message: /folded site frequency spectrum/,
  });
});

test("a draw no population can fill leaves the draw missing and the rest alone", () => {
  // A `numCalledAlleles` of 400 is every gene copy the 200 diploid
  // individuals of the panel hold, which is the largest draw the dataset
  // allows, and its missing genotypes leave no population able to call that
  // many alleles at any variant. It is not an error: the two counts of the
  // variants in a draw are 0, the three `inDraw` arrays are NaN and every bin
  // of the spectrum is 0, which is what says the question was not answered.
  const withNoDraw = ofThePanel({ stats: WITH_NO_SPECTRUM });
  const diversity = ofThePanel({
    stats: [...WITH_NO_SPECTRUM, "folded_sfs"],
    numCalledAlleles: EVERY_GENE_COPY_OF_THE_PANEL,
  });

  assert.deepEqual(diversity.numVars.inDraw, Uint32Array.of(0, 0, 0));
  assert.equal(diversity.numVarsEveryPop, PANEL_NUM_VARS);
  assert.equal(diversity.numVarsEveryPopInDraw, 0);
  const alleles = counted(diversity.numAlleles, "count of alleles");
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  const spectra = counted(diversity.foldedSfs, "folded spectrum");
  for (const [pop, name] of PANEL_POP_NAMES.entries()) {
    for (const [what, inDraw] of [
      ["alleles", alleles.inDraw],
      ["private alleles", theirPrivate.inDraw],
      ["variable variants", variable.inDraw],
    ] as const) {
      assert.ok(
        Number.isNaN(inDraw[pop] as number),
        `the ${what} of ${name} in a draw of ${EVERY_GENE_COPY_OF_THE_PANEL}`,
      );
    }
    const ours = spectra[name] as Float64Array;
    assert.equal(
      ours.length,
      EVERY_GENE_COPY_OF_THE_PANEL / 2 + 1,
      `the bins of ${name}`,
    );
    assert.ok(
      ours.every((value) => value === 0),
      `the spectrum of ${name}`,
    );
  }
  // The totals, the means, the ratio and F_IS read no draw.
  assert.deepEqual(
    alleles.total,
    counted(withNoDraw.numAlleles, "count of alleles").total,
  );
  assert.deepEqual(
    variable.ratio,
    counted(withNoDraw.variableVarsRatio, "count of variable variants").ratio,
  );
  assert.deepEqual(diversity.fis, withNoDraw.fis);
});

test("a draw larger than the dataset holds is refused and names the largest it allows", () => {
  // One allele more than every gene copy the dataset holds is a draw no
  // variant of any population could reach, so it is a user's mistake and not a
  // fact about the data. It is a different case from the draw of 400 above,
  // which this dataset's missing genotypes leave unfillable and which is no
  // error.
  assert.throws(
    () =>
      ofThePanel({
        stats: WITH_NO_SPECTRUM,
        numCalledAlleles: EVERY_GENE_COPY_OF_THE_PANEL + 1,
      }),
    {
      message:
        "`numCalledAlleles` is 401 and the largest draw this dataset allows " +
        "is 400, every gene copy of its 200 individuals at a ploidy of 2: no " +
        "population can have called more alleles than that at a variant",
    },
  );
});

test("a numCalledAlleles that is no whole number the core holds is refused by the package", () => {
  // What the package checks is that the number arrives as the number the user
  // wrote: a number of JavaScript reaches a whole number of the core as 32
  // bits with no error, so a draw of 20.5 alleles would be a draw of 20, one
  // of -1 a draw of 4294967295 and one of 4294967296 a draw of 0.
  for (const given of [20.5, -1, 4294967296]) {
    assert.throws(
      () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: given }),
      { message: /`numCalledAlleles` is a whole number of 0 or more/ },
      `numCalledAlleles: ${given}`,
    );
  }
});

test("a draw of fewer than two alleles is refused by the core", () => {
  // How small a draw may be is a rule of the core and is not written here as
  // well: a draw of one allele finds one allele whatever the population holds.
  // A Python user reads the same sentence, with `num_called_alleles` where
  // this one has the name they wrote.
  for (const given of [0, 1]) {
    assert.throws(
      () => ofThePanel({ stats: WITH_NO_SPECTRUM, numCalledAlleles: given }),
      {
        message:
          `\`numCalledAlleles\` is ${given}, and a draw shows more than one ` +
          "allele only when it is of 2 alleles at least: a draw of one " +
          "allele finds one allele whatever the population holds",
      },
      `numCalledAlleles: ${given}`,
    );
  }
});

test("the four counts of variants and the values over each of them are of the draw", () => {
  // The two variants of `SOME_VARS_SHORT_OF_THE_DRAW`, whose comment says what
  // each population calls where. The four counts of variants are four
  // different numbers here, so a value taken from the wrong one of them fails,
  // which on the panel it cannot: all four are 1200 there.
  //
  // The standardized values are the formulas of the spec on those counts, and
  // `tests/test_diversity.py` works each of them out and asserts the same
  // numbers. At the first variant `p1` called the allele 0 once and the allele
  // 1 three times of 4, so a draw of 2 shows the first with chance
  // 1 - C(3, 2) / C(4, 2) = 0.5 and the second for certain: 1.5 alleles, and
  // the same 1.5 at the second variant, where the counts are 3 and 1. `p2`
  // shows its one allele for certain at the one variant in its draw.
  //
  // The private alleles are over the one variant every population reached the
  // draw at, where `p1` holds two alleles that no draw of `p2` can show, 0.5
  // and 1, which add to 1.5, and `p2` holds one that no draw of `p1` can show,
  // 1. Over the variants of each population instead, which is the divisor of
  // the other two standardized values, `p1` would read 0.75.
  const variants = openVcf(SOME_VARS_SHORT_OF_THE_DRAW);

  const diversity = calcPopDiversity(variants, {
    pops: { p1: ["ind1", "ind2"], p2: ["ind3"] },
    stats: WITH_NO_SPECTRUM,
    numCalledAlleles: 2,
    minNumIndividuals: 0,
  });

  assert.deepEqual(diversity.pops, ["p1", "p2"]);
  assert.deepEqual(diversity.numVars.withData, Uint32Array.of(2, 2));
  assert.deepEqual(diversity.numVars.inDraw, Uint32Array.of(2, 1));
  assert.equal(diversity.numVarsEveryPop, 2);
  assert.equal(diversity.numVarsEveryPopInDraw, 1);
  const alleles = counted(diversity.numAlleles, "count of alleles");
  assert.deepEqual(alleles.total, Uint32Array.of(4, 2));
  assert.deepEqual(alleles.mean, Float64Array.of(2, 1));
  assert.deepEqual(alleles.inDraw, Float64Array.of(1.5, 1));
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  assert.deepEqual(theirPrivate.total, Uint32Array.of(3, 1));
  assert.deepEqual(theirPrivate.mean, Float64Array.of(1.5, 0.5));
  assert.deepEqual(theirPrivate.inDraw, Float64Array.of(1.5, 1));
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  assert.deepEqual(variable.total, Uint32Array.of(2, 0));
  assert.deepEqual(variable.ratio, Float64Array.of(1, 0));
  assert.deepEqual(variable.inDraw, Float64Array.of(0.5, 0));
  variants.free();
});

test("the mean private alleles are over the variants that counted for every population", () => {
  // The four variants of `SOME_VARS_MISSED`, in two populations, `p1` of
  // `ind1` and `p2` of `ind2` and `ind3`, at a threshold of one called
  // genotype. `p2` called nothing at the fourth variant, so that variant
  // counts for `p1` and not for `p2`: `p1` counts 4 variants, `p2` counts 3,
  // and 3 counted for both.
  //
  // The private alleles are over the variants that counted for every
  // population, the fourth left out, where the allele `p1` called would be
  // private only because `p2` has no data there. `p1` holds a private allele
  // at the first variant, the `0` that `p2` did not call, and at the third,
  // the `1`, which makes 2 over 3 variants, 0.6666666666666666; with the
  // variants of `p1` itself as the divisor it would be 2 over 4, 0.5. `p2`
  // holds the `1` of the first variant, 1 over 3.
  const variants = openVcf(SOME_VARS_MISSED);

  const diversity = calcPopDiversity(variants, {
    pops: { p1: ["ind1"], p2: ["ind2", "ind3"] },
    stats: ["num_alleles", "private_alleles", "variable_vars_ratio"],
    minNumIndividuals: 1,
  });

  assert.deepEqual(diversity.pops, ["p1", "p2"]);
  assert.deepEqual(diversity.numVars.withData, Uint32Array.of(4, 3));
  assert.equal(diversity.numVarsEveryPop, 3);
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  assert.deepEqual(theirPrivate.total, Uint32Array.of(2, 1));
  assert.deepEqual(
    theirPrivate.mean,
    Float64Array.of(0.6666666666666666, 0.3333333333333333),
  );
  // The alleles called and the variable variants are over the variants of
  // each population, the fourth among those of `p1`: `p1` called one allele
  // at the first, the second and the fourth and two at the third, and it is
  // the third that makes it variable.
  const alleles = counted(diversity.numAlleles, "count of alleles");
  assert.deepEqual(alleles.total, Uint32Array.of(5, 3));
  assert.deepEqual(alleles.mean, Float64Array.of(1.25, 1));
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  assert.deepEqual(variable.total, Uint32Array.of(1, 0));
  assert.deepEqual(variable.ratio, Float64Array.of(0.25, 0));
  variants.free();
});

test("a population no variant counted for has no mean, no ratio and no F_IS", () => {
  // At a threshold of 60 called genotypes the 48 individuals of `p0` can
  // never reach it, so no variant counts for that population: its totals are
  // 0, its mean, its ratio and its F_IS are NaN, and the other two
  // populations are given as they are. The variants that counted for every
  // population are 0 too, so no population has a mean private allele, `p1`
  // and `p2` included, whose own counts are above 0: that is the divisor of
  // the private alleles and not the variants of each population.
  const diversity = ofThePanel({
    stats: WITH_NO_SPECTRUM,
    minNumIndividuals: 60,
  });

  const alleles = counted(diversity.numAlleles, "count of alleles");
  const variable = counted(
    diversity.variableVarsRatio,
    "count of variable variants",
  );
  const theirPrivate = counted(
    diversity.privateAlleles,
    "count of private alleles",
  );
  const fis = diversity.fis as Float64Array;
  const p0 = diversity.pops.indexOf("p0");
  assert.equal(diversity.numVars.withData[p0], 0);
  assert.equal(alleles.total[p0], 0);
  assert.equal(variable.total[p0], 0);
  assert.ok(Number.isNaN(alleles.mean[p0] as number), "the mean of p0");
  assert.ok(Number.isNaN(variable.ratio[p0] as number), "the ratio of p0");
  assert.ok(Number.isNaN(fis[p0] as number), "the F_IS of p0");
  assert.equal(diversity.numVarsEveryPop, 0);
  for (const pop of PANEL_POP_NAMES) {
    const which = diversity.pops.indexOf(pop);
    assert.ok(
      Number.isNaN(theirPrivate.mean[which] as number),
      `the mean private alleles of ${pop}`,
    );
  }
  for (const pop of ["p1", "p2"]) {
    const which = diversity.pops.indexOf(pop);
    assert.ok((diversity.numVars.withData[which] as number) > 0, pop);
    assert.ok((alleles.total[which] as number) > 0, pop);
    assert.ok(!Number.isNaN(fis[which] as number), pop);
  }
});

test("a pass that gives no variant is refused and says where the variants went", () => {
  // The whole message is asserted, because it is the one the Python function
  // gives for the same source: the two languages say the same thing, and
  // Python writes the path of the file before it, which the bytes a
  // TypeScript user gives have not.
  const ofNoVariant = openVcf(vcfOf([]));
  try {
    assert.throws(
      () => calcPopDiversity(ofNoVariant, { stats: WITH_NO_SPECTRUM }),
      {
        message:
          "the pass gave no variant and its source holds none: a statistic " +
          "of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    ofNoVariant.free();
  }
  // A major allele frequency is at least one over the alleles of a variant,
  // so a threshold of 0 keeps none of the four variants, and the message
  // names the filter with what it was given and kept.
  const filtered = openVcf(SOME_VARS_MISSED);
  try {
    filtered.filterByMaf(0);
    assert.throws(
      () => calcPopDiversity(filtered, { stats: WITH_NO_SPECTRUM }),
      {
        message:
          "the pass gave no variant: its source gave 4 and the steps kept " +
          "none of them, the `maf` filter was given 4 and kept 0; a " +
          "statistic of a pass is calculated over the variants it gives",
      },
    );
  } finally {
    filtered.free();
  }
});

test("a population that names an individual of no pass is refused", () => {
  // `pops` is read as `calcPerVarDistribs` reads it: a name that is not an
  // individual of the pass names the population and the name.
  assert.throws(
    () =>
      ofThePanel({
        stats: WITH_NO_SPECTRUM,
        pops: { p0: ["s000", "s999"] },
      }),
    { message: /s999/ },
  );
});

test("stats names the statistics there are", () => {
  // Both messages are the core's, so a Python user reads the same sentence:
  // a name of no statistic is refused with the five there are, and a `stats`
  // that names none with what such a pass would do, read every variant of the
  // source and compute nothing of them.
  assert.throws(() => ofThePanel({ stats: ["allelic_richness"] as never }), {
    message: /`allelic_richness` is not one of the statistics/,
  });
  assert.throws(() => ofThePanel({ stats: [] }), {
    message: /`stats` names no statistic and a pass that computes none/,
  });
});

test("variants that is not a Variants is refused", () => {
  assert.throws(() => calcPopDiversity(PANEL as never), {
    message: /openVcf/,
  });
});
