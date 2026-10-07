/**
 * `calcVariantsSummary`, which gives in one pass what `calcPerVarDistribs`,
 * `calcPerIndividualStats` and `calcVarDensity` give in three.
 *
 * "The three statistics of a file in one pass" of `docs/specs/js_sources.md`
 * has the design: each of the three is what its own call gives with the
 * same options, to the bit, a statistic left out is `null`, and a call that
 * asks for none is an `Error`. The stops of a pass by `onSoFar` and the
 * options of the result so far that are refused are in `so_far.test.ts`,
 * with those of the three calls.
 *
 * The files are `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
 * individuals on chr1 and chr2, opened with the variants that passed their
 * FILTER, 475 of them, which is one block for the VCF reader; and the vars
 * file of its 500 variants written in batches of 100, which is five blocks
 * a pass, so that the result so far is given five times.
 *
 * The counts of the FILTER column, `filterColumn`, are checked against
 * `many.vcf` read with every variant: 450 with `PASS`, 25 with a dot and 25
 * with `q10`, counted with `grep -v '^#' many.vcf | cut -f7 | sort | uniq
 * -c`, so 475 passed and 25 failed. A source without the record is
 * `tests/reference/vars/of_1_1.vars`, a vars file of format 1.1.
 *
 * `assert/strict` compares two `Float64Array` byte by byte, so a mean or a
 * rate that differs in its last bit, or a NaN against a number, fails it.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type {
  PassStats,
  PerIndividualStats,
  PerVarDistribs,
  VarDensity,
  Variants,
  VariantsSummary,
} from "popnei";
import {
  calcPerIndividualStats,
  calcPerVarDistribs,
  calcVarDensity,
  calcVariantsSummary,
  init,
  openVars,
  openVcf,
  writeVars,
} from "popnei";

import { THE_POPS } from "./consumers.ts";
import { referenceVars, referenceVcf } from "./reference.ts";

await init();

/** The variants of each batch of the vars file, and of each block of a pass. */
const VARIANTS_PER_BLOCK = 100;

/** The bytes of `many.vcf`. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The bytes of the vars file of every variant of `many.vcf`, in five batches. */
const IN_FIVE_BLOCKS = await (async () => {
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  try {
    return writeVars(variants, { numVarsPerBlock: VARIANTS_PER_BLOCK }).bytes;
  } finally {
    variants.free();
  }
})();

/**
 * The options of the distributions, each away from its default so that a
 * summary that read the default in place of one of them gives other
 * numbers: two populations of 25, 5 called genotypes, 20 bins and a
 * polymorphism threshold of 0.9.
 */
const PER_VAR = {
  pops: THE_POPS,
  minNumIndividuals: 5,
  histKwargs: { numBins: 20 },
  polyThreshold: 0.9,
};

/**
 * The density in windows of 1000 base pairs over lengths past the last
 * variant of each chromosome, 30 windows of chr1 and 25 of chr2.
 */
const DENSITY = {
  windowSize: 1000,
  chromLengths: { chr1: 30000, chr2: 25000 },
};

/** The three statistics, each with its options. */
const THE_THREE = { perVar: PER_VAR, perIndividual: {}, density: DENSITY };

/** The two files, with how each is opened. */
const THE_FILES: readonly { name: string; open: () => Variants }[] = [
  {
    name: "the vars file of many.vcf in five blocks",
    open: () => openVars(IN_FIVE_BLOCKS),
  },
  { name: "many.vcf", open: () => openVcf(MANY_VCF) },
];

/** What the three calls give over `variants`, each with its options. */
function theThreeCalls(variants: Variants): {
  perVar: PerVarDistribs;
  perIndividual: PerIndividualStats;
  density: VarDensity;
} {
  return {
    perVar: calcPerVarDistribs(variants, PER_VAR),
    perIndividual: calcPerIndividualStats(variants),
    density: calcVarDensity(variants, DENSITY.windowSize, {
      chromLengths: DENSITY.chromLengths,
    }),
  };
}

/** `result` without its `passStats`, the fields that are compared whole. */
function withoutPassStats(result: { passStats: PassStats } | null): object {
  assert.ok(result !== null, "a statistic that was asked for is null");
  const { passStats: _passStats, ...theRest } = result;
  return theRest;
}

for (const file of THE_FILES) {
  test(`calcVariantsSummary over ${file.name} gives each of the three as its own call gives it, to the bit`, () => {
    const variants = file.open();
    try {
      const summary = calcVariantsSummary(variants, THE_THREE);
      const calls = theThreeCalls(variants);
      assert.deepEqual(summary.perVar, calls.perVar);
      assert.deepEqual(summary.perIndividual, calls.perIndividual);
      assert.deepEqual(summary.density, calls.density);
      assert.deepEqual(summary.passStats, calls.perVar.passStats);
      assert.deepEqual(summary.passStats, calls.perIndividual.passStats);
      assert.deepEqual(summary.passStats, calls.density.passStats);
    } finally {
      variants.free();
    }
  });

  for (const leftOut of ["perVar", "perIndividual", "density"] as const) {
    test(`calcVariantsSummary over ${file.name} without ${leftOut} gives it null and the other two unchanged`, () => {
      const variants = file.open();
      try {
        const { [leftOut]: _leftOut, ...theOtherTwo } = THE_THREE;
        const summary = calcVariantsSummary(variants, theOtherTwo);
        const calls = theThreeCalls(variants);
        assert.equal(summary[leftOut], null);
        for (const asked of Object.keys(theOtherTwo) as (keyof typeof calls)[]) {
          assert.deepEqual(summary[asked], calls[asked], asked);
        }
        assert.deepEqual(summary.passStats, calls.perVar.passStats);
      } finally {
        variants.free();
      }
    });
  }
}

test("calcVariantsSummary with soFarEvery 0 gives onSoFar the three over the variants read after each block", () => {
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    const calls: VariantsSummary[] = [];
    const result = calcVariantsSummary(variants, {
      ...THE_THREE,
      onSoFar: (soFar) => {
        calls.push(soFar);
      },
      soFarEvery: 0,
    });
    assert.deepEqual(
      calls.map((call) => call.passStats.numVars),
      [100, 200, 300, 400, 500],
    );
    for (const call of calls) {
      const numVars = call.passStats.numVars;
      // The three calls with the filter of the first n read the same first
      // variants and have one filter more, so of their counts only the
      // number of variants is the same.
      const overTheFirst = openVars(IN_FIVE_BLOCKS);
      try {
        overTheFirst.filterFirstN(numVars);
        const expected = theThreeCalls(overTheFirst);
        for (const statistic of ["perVar", "perIndividual", "density"] as const) {
          assert.deepEqual(
            withoutPassStats(call[statistic]),
            withoutPassStats(expected[statistic]),
            `the ${statistic} so far over ${numVars} variants`,
          );
          assert.equal(call[statistic]?.passStats.numVars, numVars);
          assert.equal(expected[statistic].passStats.numVars, numVars);
        }
      } finally {
        overTheFirst.free();
      }
    }
    assert.deepEqual(calls.at(-1), result);
  } finally {
    variants.free();
  }
});

test("calcVariantsSummary with none of the four is an Error that says to ask for one", () => {
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    assert.throws(() => calcVariantsSummary(variants), {
      name: "Error",
      message:
        /`calcVariantsSummary` was asked for none of its four statistics: give `perVar: \{\}`.*`perIndividual: \{\}`.*`density: \{windowSize\}`.*`filterColumn: \{\}`/,
    });
    assert.throws(
      () =>
        calcVariantsSummary(variants, {
          perVar: undefined,
          filterColumn: undefined,
        }),
      /none of its four statistics/,
    );
  } finally {
    variants.free();
  }
});

/**
 * The options of the four that are refused at the call, each with what the
 * message names: those of `perVar` and `density` are checked as the call of
 * each checks them, and `perIndividual` and `filterColumn` are empty
 * objects.
 */
const REFUSED: readonly [string, object, RegExp][] = [
  [
    "a perVar with an onSoFar of its own",
    { perVar: { onSoFar: () => {} } },
    /`onSoFar` is not an option of `calcVariantsSummary.perVar`/,
  ],
  [
    "a perVar with a stats of no statistic",
    { perVar: { stats: [] } },
    /`stats` names no statistic/,
  ],
  [
    "a perVar with a polyThreshold that is not a number",
    { perVar: { polyThreshold: "0.9" } },
    /`polyThreshold` is a number/,
  ],
  [
    "a perVar that is not an object",
    { perVar: true },
    /the options of `calcVariantsSummary.perVar` are an object/,
  ],
  [
    "a perIndividual with a key",
    { perIndividual: { pops: THE_POPS } },
    /`pops` is not an option of `calcVariantsSummary.perIndividual`, which takes no option/,
  ],
  [
    "a density with no windowSize",
    { density: {} },
    /`density.windowSize`/,
  ],
  [
    "a density with a windowSize of 0",
    { density: { windowSize: 0 } },
    /`density.windowSize`/,
  ],
  [
    "a density with chromLengths that are a Map",
    { density: { windowSize: 1000, chromLengths: new Map() } },
    /`chromLengths` is a plain object/,
  ],
  [
    "a filterColumn with a key",
    { filterColumn: { passed: true } },
    /`passed` is not an option of `calcVariantsSummary.filterColumn`, which takes no option/,
  ],
  [
    "a filterColumn that is not an object",
    { filterColumn: true },
    /the options of `calcVariantsSummary.filterColumn` are an object/,
  ],
  [
    "an option that is none of the six",
    { perVar: {}, perVars: {} },
    /`perVars` is not an option of `calcVariantsSummary`/,
  ],
];

for (const [what, options, message] of REFUSED) {
  test(`calcVariantsSummary refuses ${what} before the pass starts`, () => {
    const variants = openVars(IN_FIVE_BLOCKS);
    try {
      let told = 0;
      variants.onProgress(() => {
        told += 1;
      });
      assert.throws(
        () => calcVariantsSummary(variants, options as never),
        message,
      );
      // A pass that started would have told the page of its first read.
      assert.equal(told, 0);
    } finally {
      variants.free();
    }
  });
}

/** `many.vcf` opened with every variant, those that failed their FILTER too. */
function everyVariantOfMany(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/**
 * The two sources of every variant of `many.vcf`, each of which gives the
 * counts of the FILTER column: the VCF itself and the vars file written
 * from it.
 */
const EVERY_VARIANT: readonly { name: string; open: () => Variants }[] = [
  { name: "many.vcf with every variant", open: everyVariantOfMany },
  {
    name: "the vars file of every variant of many.vcf",
    open: () => openVars(IN_FIVE_BLOCKS),
  },
];

for (const file of EVERY_VARIANT) {
  test(`calcVariantsSummary with filterColumn alone over ${file.name} counts 475 passed and 25 failed in a pass of 500`, () => {
    const variants = file.open();
    try {
      const summary = calcVariantsSummary(variants, { filterColumn: {} });
      assert.deepEqual(summary.filterColumn, { passed: 475, failed: 25 });
      assert.equal(summary.passStats.numVars, 500);
      assert.equal(summary.perVar, null);
      assert.equal(summary.perIndividual, null);
      assert.equal(summary.density, null);
    } finally {
      variants.free();
    }
  });

  test(`calcVariantsSummary with filterColumn beside the three over ${file.name} counts 475 and 25 and leaves the three as they are without it`, () => {
    const variants = file.open();
    try {
      const withTheCounts = calcVariantsSummary(variants, {
        ...THE_THREE,
        filterColumn: {},
      });
      const without = calcVariantsSummary(variants, THE_THREE);
      assert.deepEqual(withTheCounts.filterColumn, { passed: 475, failed: 25 });
      assert.equal(without.filterColumn, null);
      assert.deepEqual(withTheCounts.perVar, without.perVar);
      assert.deepEqual(withTheCounts.perIndividual, without.perIndividual);
      assert.deepEqual(withTheCounts.density, without.density);
      assert.deepEqual(withTheCounts.passStats, without.passStats);
      assert.equal(withTheCounts.passStats.numVars, 500);
    } finally {
      variants.free();
    }
  });

  test(`calcVariantsSummary with filterColumn over ${file.name} after filterPassed counts 475 passed and none failed`, () => {
    const variants = file.open();
    try {
      variants.filterPassed();
      const summary = calcVariantsSummary(variants, { filterColumn: {} });
      assert.deepEqual(summary.filterColumn, { passed: 475, failed: 0 });
      assert.equal(summary.passStats.numVars, 475);
      // The variants the filter took out are in its own counts.
      assert.deepEqual(summary.passStats.filtering.passed, {
        varsProcessed: 500,
        varsKept: 475,
      });
    } finally {
      variants.free();
    }
  });
}

test("calcVariantsSummary with filterColumn and soFarEvery 0 gives onSoFar the counts over the variants read after each block", () => {
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    const calls: VariantsSummary[] = [];
    const result = calcVariantsSummary(variants, {
      filterColumn: {},
      onSoFar: (soFar) => {
        calls.push(soFar);
      },
      soFarEvery: 0,
    });
    assert.deepEqual(
      calls.map((call) => call.passStats.numVars),
      [100, 200, 300, 400, 500],
    );
    for (const call of calls) {
      const numVars = call.passStats.numVars;
      const overTheFirst = openVars(IN_FIVE_BLOCKS);
      try {
        overTheFirst.filterFirstN(numVars);
        const expected = calcVariantsSummary(overTheFirst, {
          filterColumn: {},
        });
        assert.ok(call.filterColumn !== null, `no counts over ${numVars}`);
        assert.deepEqual(
          call.filterColumn,
          expected.filterColumn,
          `the counts so far over ${numVars} variants`,
        );
        assert.equal(
          call.filterColumn.passed + call.filterColumn.failed,
          numVars,
        );
      } finally {
        overTheFirst.free();
      }
    }
    assert.deepEqual(calls.at(-1), result);
    assert.deepEqual(result.filterColumn, { passed: 475, failed: 25 });
  } finally {
    variants.free();
  }
});

/** The message of a source that holds no record of the FILTER column. */
const NOT_RECORDED =
  /the variants hold no record of whether they passed their FILTER, so the summary cannot count how many passed and how many failed/;

test("calcVariantsSummary with filterColumn over a vars file of 1.1, which has no record of whether its variants passed, is an Error", async () => {
  const variants = openVars(await referenceVars("of_1_1.vars"));
  try {
    assert.equal(variants.keepsPassed, false);
    assert.throws(() => calcVariantsSummary(variants, { filterColumn: {} }), {
      name: "Error",
      message: NOT_RECORDED,
    });
    assert.throws(
      () =>
        calcVariantsSummary(variants, { perIndividual: {}, filterColumn: {} }),
      NOT_RECORDED,
    );
  } finally {
    variants.free();
  }
});

test("calcVariantsSummary with filterColumn over a vars file of 1.1 gives the errors of perVar and density before that of the record", async () => {
  const variants = openVars(await referenceVars("of_1_1.vars"));
  try {
    assert.throws(
      () =>
        calcVariantsSummary(variants, {
          perVar: { polyThreshold: 2 },
          density: { windowSize: 1000 },
          filterColumn: {},
        }),
      /`polyThreshold` is 2, and a threshold is a number from 0 to 1/,
    );
    // Windows of 1 base pair over 20000000 base pairs are more than the
    // 10000000 the core gives at most, which only the core refuses.
    assert.throws(
      () =>
        calcVariantsSummary(variants, {
          density: { windowSize: 1, chromLengths: { chr1: 20_000_000 } },
          filterColumn: {},
        }),
      /the density of the variants in windows of 1 base pairs has 20000000 windows at least/,
    );
  } finally {
    variants.free();
  }
});
