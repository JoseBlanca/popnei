/**
 * The result so far of the three calculations whose results add up over the
 * blocks of a pass: `calcPerVarDistribs`, `calcPerIndividualStats` and
 * `calcVarDensity`, and of `calcVariantsSummary`, which gives the three in
 * one pass, given to the function `onSoFar` while their pass runs.
 *
 * "The result so far" of `docs/specs/js_sources.md` has the design. The
 * function is called after a block, the last one too, when `soFarEvery`
 * seconds have gone by since the pass started or since the last call, with
 * the result the consumer would return over the variants of the blocks read
 * so far. A value it throws ends the pass as one that the function of
 * `onProgress` stopped: the consumer throws that value, and the page is told
 * nothing more of how far the pass read.
 *
 * The file is `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50
 * individuals on chr1, from 1000 to 10213, and chr2, from 10250 to 19463,
 * opened with every variant and written as a vars file in batches of 100
 * variants, so that each pass is five blocks of 100. `many.vcf` itself is
 * one block for the VCF reader, which cuts its blocks by the genotypes they
 * hold, and would give one call.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Progress, Variants } from "popnei";
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

import { referenceVcf } from "./reference.ts";

await init();

/** How many variants `many.vcf` holds, every one of them in the vars file. */
const VARIANTS_OF_MANY_VCF = 500;

/** The variants of each batch of the vars file, and of each block of a pass. */
const VARIANTS_PER_BLOCK = 100;

/** The bytes of the vars file of `many.vcf`, in five batches of 100. */
const IN_FIVE_BLOCKS = await (async () => {
  const variants = openVcf(await referenceVcf("many.vcf"), {
    onlyPassed: false,
  });
  try {
    return writeVars(variants, { numVarsPerBlock: VARIANTS_PER_BLOCK }).bytes;
  } finally {
    variants.free();
  }
})();

/** The width of the windows of the density, in base pairs. */
const WINDOW_SIZE = 1000;

/**
 * Lengths for the two chromosomes of `many.vcf`, past their last variants,
 * which lay out every window from the first call: 30 of chr1 and 25 of chr2.
 */
const CHROM_LENGTHS = { chr1: 30000, chr2: 25000 };

/** What every result of the four carries. */
interface WithPassStats {
  passStats: PassStats;
}

/** The two options of the result so far, as every one of the four takes them. */
interface SoFarOptions {
  onSoFar?: unknown;
  soFarEvery?: unknown;
}

/** One of the four calculations, with a call of it over `variants`. */
interface TheCalculation {
  /** What the names of the tests call it. */
  name: string;
  /** The call, with the two options of the result so far. */
  run: (variants: Variants, options: SoFarOptions) => WithPassStats;
}

/**
 * The four calculations, the density twice: with lengths, whose windows are
 * all there from the first call, and with `chromLengths: {}`, which takes no
 * length and has the windows grow up to the last variant read so far.
 * `many.vcf` has no `##contig` length, so its own lengths would give the
 * second case again. `calcVariantsSummary` gives the three, the density with
 * lengths; `variants_summary.test.ts` compares its result so far with those
 * of the three calls.
 *
 * The options cross as `never`: the tests of the arguments give values that
 * the types of the package refuse, which is what a user of JavaScript does.
 */
const THE_CALCULATIONS: readonly TheCalculation[] = [
  {
    name: "calcPerVarDistribs",
    run: (variants, options) => calcPerVarDistribs(variants, options as never),
  },
  {
    name: "calcPerIndividualStats",
    run: (variants, options) =>
      calcPerIndividualStats(variants, options as never),
  },
  {
    name: "calcVarDensity with chromLengths",
    run: (variants, options) =>
      calcVarDensity(variants, WINDOW_SIZE, {
        chromLengths: CHROM_LENGTHS,
        ...(options as object),
      }),
  },
  {
    name: "calcVarDensity with chromLengths {}",
    run: (variants, options) =>
      calcVarDensity(variants, WINDOW_SIZE, {
        chromLengths: {},
        ...(options as object),
      }),
  },
  {
    name: "calcVariantsSummary",
    run: (variants, options) =>
      calcVariantsSummary(variants, {
        perVar: {},
        perIndividual: {},
        density: { windowSize: WINDOW_SIZE, chromLengths: CHROM_LENGTHS },
        ...(options as object),
      }),
  },
];

/**
 * `result` without its `passStats`, the fields that are compared whole, and
 * without those of the three results a `calcVariantsSummary` holds.
 */
function withoutPassStats(result: WithPassStats): object {
  const { passStats: _passStats, ...theRest } = result;
  return Object.fromEntries(
    Object.entries(theRest).map(([key, value]) => [
      key,
      typeof value === "object" && value !== null && "passStats" in value
        ? withoutPassStats(value as WithPassStats)
        : value,
    ]),
  );
}

/** What `calculation` gives over the first `numVars` variants of the file. */
function overTheFirst(
  calculation: TheCalculation,
  numVars: number,
): WithPassStats {
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    variants.filterFirstN(numVars);
    return calculation.run(variants, {});
  } finally {
    variants.free();
  }
}

/**
 * What `run` threw, and an assertion that fails when it threw nothing.
 *
 * What these tests stop a pass with is not an `Error`, and what they assert
 * of it is that it is one value and not another that reads the same.
 */
function whatWasThrownBy(run: () => void): unknown {
  try {
    run();
  } catch (thrown: unknown) {
    return thrown;
  }
  return assert.fail("the run gave its result instead of being stopped");
}

/**
 * What the function of these tests throws, which is not an `Error`: an
 * application stops a pass with a value of its own and tells it from a file
 * that could not be read by what it is.
 */
const THE_CANCEL = { whyThePassEnded: "the user pressed the button" };

for (const calculation of THE_CALCULATIONS) {
  test(`${calculation.name} with soFarEvery 0 calls onSoFar after each block with the result over the variants read`, () => {
    const variants = openVars(IN_FIVE_BLOCKS);
    try {
      const calls: WithPassStats[] = [];
      const result = calculation.run(variants, {
        onSoFar: (soFar: WithPassStats) => {
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
        // The pass with the filter of the first n reads the same first
        // variants and has one filter more, so of its counts only the
        // number of variants is the same.
        const expected = overTheFirst(calculation, numVars);
        assert.deepEqual(
          withoutPassStats(call),
          withoutPassStats(expected),
          `the result so far over ${numVars} variants`,
        );
        assert.equal(expected.passStats.numVars, numVars);
      }
      assert.deepEqual(calls.at(-1), result);
      assert.equal(result.passStats.numVars, VARIANTS_OF_MANY_VCF);
    } finally {
      variants.free();
    }
  });
}

test("calcVarDensity with chromLengths gives every window at the first onSoFar call", () => {
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    const numWindows: number[] = [];
    calcVarDensity(variants, WINDOW_SIZE, {
      chromLengths: CHROM_LENGTHS,
      onSoFar: (soFar) => {
        numWindows.push(soFar.chroms.length);
      },
      soFarEvery: 0,
    });
    assert.deepEqual(numWindows, [55, 55, 55, 55, 55]);
  } finally {
    variants.free();
  }
});

test("calcVarDensity with chromLengths {} gives onSoFar the windows up to the last variant read", () => {
  // The first 100 variants are on chr1 alone, and so are the first 250.
  const variants = openVars(IN_FIVE_BLOCKS);
  try {
    const chromsCalledWith: string[][] = [];
    calcVarDensity(variants, WINDOW_SIZE, {
      chromLengths: {},
      onSoFar: (soFar) => {
        chromsCalledWith.push([...new Set(soFar.chroms)]);
      },
      soFarEvery: 0,
    });
    assert.deepEqual(chromsCalledWith.at(0), ["chr1"]);
    assert.deepEqual(chromsCalledWith.at(-1), ["chr1", "chr2"]);
  } finally {
    variants.free();
  }
});

for (const calculation of THE_CALCULATIONS) {
  test(`${calculation.name} with soFarEvery 3600 never calls onSoFar over a file of seconds`, () => {
    const variants = openVars(IN_FIVE_BLOCKS);
    try {
      let calls = 0;
      const result = calculation.run(variants, {
        onSoFar: () => {
          calls += 1;
        },
        soFarEvery: 3600,
      });
      assert.equal(calls, 0);
      assert.equal(result.passStats.numVars, VARIANTS_OF_MANY_VCF);
    } finally {
      variants.free();
    }
  });
}

for (const calculation of THE_CALCULATIONS) {
  test(`${calculation.name} throws the value onSoFar threw at its second call, and the next call runs whole`, () => {
    const variants = openVars(IN_FIVE_BLOCKS);
    try {
      // The page is told of the first read of the pass and of nothing after
      // the stop: the vars file is smaller than one range, so the only other
      // call would be the one at the end of the run, which says how far the
      // pass read and which a pass that was stopped is not given.
      const told: Progress[] = [];
      variants.onProgress((progress) => {
        told.push(progress);
      });
      const calls: number[] = [];
      const thrown = whatWasThrownBy(() => {
        calculation.run(variants, {
          onSoFar: (soFar: WithPassStats) => {
            calls.push(soFar.passStats.numVars);
            if (calls.length === 2) {
              throw THE_CANCEL;
            }
          },
          soFarEvery: 0,
        });
      });
      assert.equal(thrown, THE_CANCEL);
      assert.deepEqual(calls, [100, 200]);
      assert.deepEqual(
        told.map((progress) => progress.bytesRead),
        [0],
      );
      // The same `Variants` runs the same calculation again from the start
      // of the file and to its end.
      told.length = 0;
      const numVarsCalledWith: number[] = [];
      const result = calculation.run(variants, {
        onSoFar: (soFar: WithPassStats) => {
          numVarsCalledWith.push(soFar.passStats.numVars);
        },
        soFarEvery: 0,
      });
      assert.deepEqual(numVarsCalledWith, [100, 200, 300, 400, 500]);
      assert.equal(result.passStats.numVars, VARIANTS_OF_MANY_VCF);
      assert.equal(told.length, 2);
    } finally {
      variants.free();
    }
  });
}

for (const calculation of THE_CALCULATIONS) {
  test(`free from inside onSoFar while ${calculation.name} reads is refused`, () => {
    const variants = openVars(IN_FIVE_BLOCKS);
    try {
      const thrown = whatWasThrownBy(() => {
        calculation.run(variants, {
          onSoFar: () => {
            variants.free();
          },
          soFarEvery: 0,
        });
      });
      assert.ok(thrown instanceof Error, `it threw ${String(thrown)}`);
      assert.match(thrown.message, /a run is reading these variants/);
      let numVars = 0;
      for (const block of variants.iterBlocks()) {
        numVars += block.numVars;
      }
      assert.equal(numVars, VARIANTS_OF_MANY_VCF);
    } finally {
      variants.free();
    }
  });
}

/**
 * The options each of the four refuses at the call, with what the message
 * names.
 */
const REFUSED: readonly [string, SoFarOptions, RegExp][] = [
  ["a soFarEvery with no onSoFar", { soFarEvery: 1 }, /`soFarEvery`.*`onSoFar`/],
  ["an onSoFar of 3", { onSoFar: 3 }, /`onSoFar` is a function/],
  ["a soFarEvery of -1", { onSoFar: () => {}, soFarEvery: -1 }, /`soFarEvery`/],
  [
    "a soFarEvery of NaN",
    { onSoFar: () => {}, soFarEvery: Number.NaN },
    /`soFarEvery`/,
  ],
  [
    "a soFarEvery of Infinity",
    { onSoFar: () => {}, soFarEvery: Number.POSITIVE_INFINITY },
    /`soFarEvery`/,
  ],
  [
    "a soFarEvery of the string 2",
    { onSoFar: () => {}, soFarEvery: "2" },
    /`soFarEvery`/,
  ],
];

for (const calculation of THE_CALCULATIONS) {
  for (const [what, options, message] of REFUSED) {
    test(`${calculation.name} refuses ${what} for onSoFar before the pass starts`, () => {
      const variants = openVars(IN_FIVE_BLOCKS);
      try {
        // A pass that started would have told the page of its first read.
        const told: Progress[] = [];
        variants.onProgress((progress) => {
          told.push(progress);
        });
        assert.throws(() => calculation.run(variants, options), message);
        assert.equal(told.length, 0);
      } finally {
        variants.free();
      }
    });
  }
}
