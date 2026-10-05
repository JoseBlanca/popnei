/**
 * The filter of the first n variants from TypeScript, `filterFirstN`: which
 * variants it keeps, the counts of a pass it ended, what it refuses, and
 * that the counts of every consumer say whether it ended their pass.
 *
 * `docs/specs/filters.md` has the filter, in "The filter that keeps the
 * first n variants". The file is `many.vcf` of `docs/specs/io_vcf.md`, 500
 * variants of 50 diploid individuals, read with every variant given, those
 * that failed their FILTER too. The positions are those bcftools 1.24 gives,
 * `bcftools view -H many.vcf | head -n 10` and, with the MAF filter of 0.8
 * before, `bcftools view -H -Q 0.8:major many.vcf | head -n 10`, which the
 * spec has from 5 October 2026.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { ConsumerName, PassStats, Variants } from "popnei";
import {
  calcGwas,
  calcKinship,
  calcLdAndDistPerPop,
  calcPairwiseKosmanDists,
  calcPerIndividualStats,
  calcPerVarDistribs,
  calcPopDists,
  calcPopDiversity,
  calcRogersHuffR2Matrix,
  calcVarDensity,
  doPcaFromVariants,
  doPcoaFromVariants,
  init,
  openVcf,
  writeVars,
  writeVcf,
} from "popnei";

import { THE_POPS, THE_TRAIT, theConsumersTheCrateNames } from "./consumers.ts";
import { referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The variants of `many.vcf`, read with every variant given. */
const MANY_NUM_VARS = 500;

/** The first ten variants of `many.vcf`, by position, all of chr1. */
const THE_FIRST_TEN = [
  1000, 1037, 1074, 1111, 1148, 1185, 1222, 1259, 1296, 1333,
];

/** The first ten that the MAF filter of 0.8 keeps, by position. */
const THE_FIRST_TEN_AFTER_THE_MAF = [
  1037, 1074, 1111, 1148, 1222, 1259, 1296, 1333, 1370, 1407,
];

/** The 500 variants of `many.vcf`, the ones that failed their FILTER among
 * them. */
function many(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/** The position of every variant one whole pass over `variants` gives, in
 * blocks of `numVarsPerBlock`, and the counts of that pass. */
function keptBy(
  variants: Variants,
  numVarsPerBlock?: number,
): { positions: number[]; passStats: PassStats } {
  const blocks = variants.iterBlocks({ fields: ["pos"], numVarsPerBlock });
  const positions: number[] = [];
  for (const block of blocks) {
    if (block.pos === null) {
      throw new Error("the pass was asked for the positions and gave none");
    }
    positions.push(...block.pos);
  }
  return { positions, passStats: blocks.passStats };
}

test("filterFirstN of 10 keeps the first ten of bcftools and ends the pass", () => {
  for (const numVarsPerBlock of [7, undefined]) {
    const variants = many();
    assert.equal(variants.filterFirstN(10), undefined);

    const { positions, passStats } = keptBy(variants, numVarsPerBlock);

    assert.deepEqual(positions, THE_FIRST_TEN);
    assert.equal(passStats.numVars, 10);
    assert.equal(passStats.stoppedEarly, true);
    assert.deepEqual(Object.keys(passStats.filtering), ["first_n"]);
    assert.equal(passStats.filtering["first_n"]?.varsKept, 10);
    variants.free();
  }
});

test("filterFirstN in blocks of 7 is given the two blocks that hold the ten", () => {
  // The filter takes whole blocks and keeps of the second only the three it
  // needs, so its counts are of the 14 variants of two blocks.
  const variants = many();
  variants.filterFirstN(10);

  const { passStats } = keptBy(variants, 7);

  assert.deepEqual(passStats, {
    numVars: 10,
    filtering: { first_n: { varsProcessed: 14, varsKept: 10 } },
    stoppedEarly: true,
  });
  variants.free();
});

test("filterFirstN after the MAF filter of 0.8 keeps the first ten it keeps", () => {
  for (const numVarsPerBlock of [7, undefined]) {
    const variants = many();
    variants.filterByMaf(0.8);
    variants.filterFirstN(10);

    const { positions, passStats } = keptBy(variants, numVarsPerBlock);

    assert.deepEqual(positions, THE_FIRST_TEN_AFTER_THE_MAF);
    assert.equal(passStats.stoppedEarly, true);
    if (numVarsPerBlock === 7) {
      // The MAF filter keeps 12 of the first 14 variants, two blocks of 7,
      // and the filter of the first n keeps 10 of those 12.
      assert.deepEqual(passStats.filtering, {
        maf: { varsProcessed: 14, varsKept: 12 },
        first_n: { varsProcessed: 12, varsKept: 10 },
      });
    }
    variants.free();
  }
});

test("filterFirstN of 1000 over 500 variants keeps them all and did not stop early", () => {
  const variants = many();
  variants.filterFirstN(1000);

  const { positions, passStats } = keptBy(variants);

  assert.equal(positions.length, MANY_NUM_VARS);
  assert.deepEqual(passStats, {
    numVars: MANY_NUM_VARS,
    filtering: {
      first_n: { varsProcessed: MANY_NUM_VARS, varsKept: MANY_NUM_VARS },
    },
    stoppedEarly: false,
  });
  variants.free();
});

test("filterFirstN adds the step first_n with numVars", () => {
  const variants = many();
  variants.filterByMaf(0.8);
  variants.filterFirstN(1000);
  assert.deepEqual(variants.steps, [
    { kind: "maf", args: { maxAllowedMaf: 0.8 } },
    { kind: "first_n", args: { numVars: 1000 } },
  ]);
  variants.free();
});

test("filterFirstN of 0 is an Error that asks for 1 or more", () => {
  const variants = many();
  assert.throws(() => variants.filterFirstN(0), {
    name: "Error",
    message: /asked for 0 variants.*ask for 1 or more/,
  });
  assert.deepEqual(variants.steps, []);
  variants.free();
});

for (const [given, written] of [
  [-1, "-1"],
  [1.5, "1.5"],
  [Number.NaN, "NaN"],
  [2 ** 53, "9007199254740992"],
  [true, "the boolean true"],
  ["10", "the string `10`"],
  [undefined, "undefined"],
] as [unknown, string][]) {
  test(`filterFirstN of ${written} is an Error that names numVars`, () => {
    // The generated code would turn each of them into a float64 with no
    // error, `true` into 1 and `undefined` into NaN, so the package refuses
    // them before the call.
    const variants = many();
    assert.throws(() => variants.filterFirstN(given as number), {
      name: "Error",
      message: new RegExp(`\`numVars\` is a whole number of variants.*${written}`),
    });
    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

/** Each method that adds a filter of the variants, under its kind, with a
 * call of it. */
const FILTERS_OF_THE_VARIANTS: [string, (variants: Variants) => void][] = [
  ["missing_data", (variants) => variants.filterByMissingData(0.1)],
  ["maf", (variants) => variants.filterByMaf(0.8)],
  ["obs_het", (variants) => variants.filterByObsHet(0.5)],
  ["ld", (variants) => variants.filterByLd(0.3, 10000)],
  [
    "regions",
    (variants) => variants.filterByRegions(new TextEncoder().encode("chr1\t0\t2000\n")),
  ],
  [
    "excluded_regions",
    (variants) =>
      variants.filterByRegions(new TextEncoder().encode("chr1\t0\t2000\n"), {
        exclude: true,
      }),
  ],
];

for (const [kind, filter] of FILTERS_OF_THE_VARIANTS) {
  test(`filterFirstN refuses a filter by ${kind} added after it`, () => {
    const variants = many();
    variants.filterFirstN(10);
    assert.throws(() => filter(variants), {
      name: "Error",
      message: new RegExp(
        `filtered by first_n already, and a filter by ${kind} after it`,
      ),
    });
    assert.deepEqual(variants.steps, [
      { kind: "first_n", args: { numVars: 10 } },
    ]);
    variants.free();
  });
}

test("filterFirstN then a second MAF filter is refused as a second of its kind", () => {
  // The MAF filter of 0.8 breaks both rules, and gets the refusal of a second
  // filter of its kind, with both thresholds, as in Python.
  const variants = many();
  variants.filterByMaf(0.9);
  variants.filterFirstN(10);
  assert.throws(() => variants.filterByMaf(0.8), {
    name: "Error",
    message:
      /filtered by maf already, with a threshold of 0\.9, and a second filter of that kind, whose threshold is 0\.8/,
  });
  assert.deepEqual(
    variants.steps.map((step) => step.kind),
    ["maf", "first_n"],
  );
  variants.free();
});

test("filterFirstN accepts the filter of individuals after it", () => {
  const variants = many();
  variants.filterFirstN(10);
  variants.filterIndividuals(["ind03", "ind01"]);

  const { positions, passStats } = keptBy(variants);

  assert.deepEqual(variants.individuals, ["ind03", "ind01"]);
  assert.deepEqual(positions, THE_FIRST_TEN);
  assert.equal(passStats.stoppedEarly, true);
  variants.free();
});

test("filterFirstN refuses a second filter of the first n", () => {
  const variants = many();
  variants.filterFirstN(10);
  assert.throws(() => variants.filterFirstN(20), {
    name: "Error",
    message: /filtered by first_n already, and a second filter of that kind/,
  });
  variants.free();
});

/**
 * The counts of each consumer, called over `variants`.
 *
 * The calls are those of `consumers.ts`, which gives none of the results
 * back; this one is held to the consumers of the crate by the test below,
 * as that list is.
 */
const THE_COUNTS_OF_EACH_CONSUMER: Record<
  ConsumerName,
  (variants: Variants) => PassStats | undefined
> = {
  calcPerVarDistribs: (variants) => calcPerVarDistribs(variants).passStats,
  calcPerIndividualStats: (variants) =>
    calcPerIndividualStats(variants).passStats,
  calcVarDensity: (variants) => calcVarDensity(variants, 1000).passStats,
  calcPairwiseKosmanDists: (variants) =>
    calcPairwiseKosmanDists(variants).passStats,
  calcPopDists: (variants) =>
    calcPopDists(variants, THE_POPS, {
      measures: ["fst"],
      jackknifeGroup: null,
      minNumIndividuals: 1,
    }).passStats,
  calcPopDiversity: (variants) =>
    calcPopDiversity(variants, { pops: THE_POPS, minNumIndividuals: 1 })
      .passStats,
  calcRogersHuffR2Matrix: (variants) =>
    calcRogersHuffR2Matrix(variants).passStats,
  calcLdAndDistPerPop: (variants) =>
    calcLdAndDistPerPop(variants, { pops: THE_POPS }).passStats,
  calcKinship: (variants) =>
    calcKinship(variants, { transformToBiallelic: true }).passStats,
  doPcaFromVariants: (variants) =>
    doPcaFromVariants(variants, {
      numPrinComps: 10,
      transformToBiallelic: true,
    }).passStats,
  doPcoaFromVariants: (variants) =>
    doPcoaFromVariants(variants, { correctByLingoes: true }).passStats,
  calcGwas: (variants) =>
    calcGwas(variants, {
      phenotype: THE_TRAIT,
      trait: "continuous",
      transformToBiallelic: true,
    }).passStats,
  writeVars: (variants) => writeVars(variants).passStats,
  writeVcf: (variants) => writeVcf(variants).passStats,
  iterBlocks: (variants) => {
    const blocks = variants.iterBlocks();
    for (const _block of blocks) {
      // Every block is read, so the pass ends where the filter ends it.
    }
    return blocks.passStats;
  },
};

test("filterFirstN is told by stoppedEarly in the counts of every consumer", () => {
  assert.deepEqual(
    Object.keys(THE_COUNTS_OF_EACH_CONSUMER).sort(),
    theConsumersTheCrateNames().sort(),
  );
  for (const [name, countsOf] of Object.entries(THE_COUNTS_OF_EACH_CONSUMER)) {
    const stopped = many();
    stopped.filterFirstN(100);
    const ofTheStopped = countsOf(stopped);
    assert.equal(ofTheStopped?.stoppedEarly, true, name);
    assert.equal(ofTheStopped?.filtering["first_n"]?.varsKept, 100, name);
    stopped.free();

    const whole = many();
    whole.filterFirstN(1000);
    const ofTheWhole = countsOf(whole);
    assert.equal(ofTheWhole?.stoppedEarly, false, name);
    assert.equal(ofTheWhole?.filtering["first_n"]?.varsKept, MANY_NUM_VARS, name);
    whole.free();
  }
});
