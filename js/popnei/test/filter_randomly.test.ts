/**
 * The filter that keeps variants at random from TypeScript,
 * `filterRandomly`: which variants it keeps for a keep rate and a seed,
 * that a second pass keeps the same ones, its counts and step, and what it
 * refuses.
 *
 * `docs/specs/filters.md` has the filter, in "The filter that keeps
 * variants at random". The file is `many.vcf` of `docs/specs/io_vcf.md`,
 * 500 variants of 50 diploid individuals, read with every variant given,
 * those that failed their FILTER too, as the numbers of the spec are of the
 * 500. The positions are those of the table of the spec, which
 * `tests/reference/filters/random_draws.py` gave on 5 October 2026 from the
 * draws of SplitMix64 that Java's `SplittableRandom` gives.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Variants } from "popnei";
import { init, openVcf } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The variants of `many.vcf`, read with every variant given. */
const MANY_NUM_VARS = 500;

/** How many variants a keep rate of 0.1 and a seed of 42 keep of the 500. */
const KEPT_AT_0_1_SEED_42 = 45;

/** The first five of them, by position, all of chr1. */
const FIRST_FIVE_AT_0_1_SEED_42 = [1148, 1666, 1777, 1888, 2332];

/** How many a keep rate of 0.1 and a seed of 7 keep, and the first five. */
const KEPT_AT_0_1_SEED_7 = 49;
const FIRST_FIVE_AT_0_1_SEED_7 = [1037, 1962, 2147, 2332, 2591];

/** The first ten that a keep rate of 0.5 and a seed of 42 keep. */
const FIRST_TEN_AT_0_5_SEED_42 = [
  1037, 1074, 1111, 1148, 1222, 1296, 1370, 1407, 1555, 1592,
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

test("filterRandomly at 0.1 and a seed of 42 keeps the 45 of the spec in blocks of 7 and of the default size", () => {
  for (const numVarsPerBlock of [7, undefined]) {
    const variants = many();
    assert.equal(variants.filterRandomly(0.1, { seed: 42 }), undefined);

    const { positions, passStats } = keptBy(variants, numVarsPerBlock);

    assert.equal(positions.length, KEPT_AT_0_1_SEED_42);
    assert.deepEqual(positions.slice(0, 5), FIRST_FIVE_AT_0_1_SEED_42);
    assert.deepEqual(passStats, {
      numVars: KEPT_AT_0_1_SEED_42,
      filtering: {
        random: { varsProcessed: MANY_NUM_VARS, varsKept: KEPT_AT_0_1_SEED_42 },
      },
      stoppedEarly: false,
    });
    variants.free();
  }
});

test("filterRandomly keeps the same variants in a second pass", () => {
  const variants = many();
  variants.filterRandomly(0.1, { seed: 42 });

  const first = keptBy(variants);
  const second = keptBy(variants);

  assert.equal(first.positions.length, KEPT_AT_0_1_SEED_42);
  assert.deepEqual(second.positions, first.positions);
  variants.free();
});

test("filterRandomly with a seed of 7 keeps the 49 of the spec", () => {
  const variants = many();
  variants.filterRandomly(0.1, { seed: 7 });

  const { positions } = keptBy(variants);

  assert.equal(positions.length, KEPT_AT_0_1_SEED_7);
  assert.deepEqual(positions.slice(0, 5), FIRST_FIVE_AT_0_1_SEED_7);
  variants.free();
});

test("filterRandomly with no seed keeps what a seed of 42 keeps and says the seed in its step", () => {
  const withTheDefault = many();
  withTheDefault.filterRandomly(0.1);
  const withSeed42 = many();
  withSeed42.filterRandomly(0.1, { seed: 42 });

  assert.deepEqual(keptBy(withTheDefault).positions, keptBy(withSeed42).positions);
  assert.deepEqual(withTheDefault.steps, [
    { kind: "random", args: { keepRate: 0.1, seed: 42 } },
  ]);
  withTheDefault.free();
  withSeed42.free();
});

test("filterRandomly with a seed of 2 ** 53 - 1 is accepted and read back whole", () => {
  // The largest whole number a number of JavaScript counts to one by one,
  // which crosses as a float64 and comes back from the core unchanged.
  const variants = many();
  variants.filterRandomly(0.1, { seed: 2 ** 53 - 1 });
  assert.deepEqual(variants.steps, [
    { kind: "random", args: { keepRate: 0.1, seed: 9007199254740991 } },
  ]);
  variants.free();
});

test("filterRandomly at 0.5 and then filterFirstN of 10 keep the ten of the spec", () => {
  const variants = many();
  variants.filterRandomly(0.5);
  variants.filterFirstN(10);

  const { positions, passStats } = keptBy(variants);

  assert.deepEqual(positions, FIRST_TEN_AT_0_5_SEED_42);
  assert.equal(passStats.stoppedEarly, true);
  variants.free();
});

for (const [given, written] of [
  [1.5, "1.5"],
  [-0.1, "-0.1"],
  [Number.NaN, "NaN"],
  [2, "2"],
  [Number.POSITIVE_INFINITY, "Infinity"],
] as [number, string][]) {
  test(`filterRandomly with a keep rate of ${written} is an Error that names keepRate`, () => {
    // The core refuses it and names the argument as Python writes it, which
    // the binding writes as TypeScript does, and the number as JavaScript
    // writes it: 2 and not the 2.0 of Rust, Infinity and not inf.
    const variants = many();
    assert.throws(() => variants.filterRandomly(given), (error: unknown) => {
      assert.ok(error instanceof Error);
      assert.ok(
        error.message.includes(`\`keepRate\` is ${written}, `),
        error.message,
      );
      assert.match(error.message, /a number from 0 to 1, both included/);
      assert.doesNotMatch(error.message, /keep_rate/);
      return true;
    });
    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

test("filterRandomly with a keep rate that is no number is an Error that names keepRate", () => {
  // The generated code would hand the core 0.1 for the string "0.1" and NaN
  // for `undefined`, with no error.
  const variants = many();
  assert.throws(() => variants.filterRandomly("0.1" as unknown as number), {
    name: "Error",
    message: /`keepRate` is a number, and the string `0\.1` was given/,
  });
  assert.deepEqual(variants.steps, []);
  variants.free();
});

for (const [given, written] of [
  [-1, "-1"],
  [1.5, "1.5"],
  [2 ** 53, "9007199254740992"],
  [Number.NaN, "NaN"],
  [null, "null"],
  ["7", "the string `7`"],
] as [unknown, string][]) {
  test(`filterRandomly with a seed of ${written} is an Error that names seed`, () => {
    // The package refuses them before the call, with the argument and what
    // was given; the binding crate would refuse those that cross as a
    // float64 as a defect of popnei.
    const variants = many();
    assert.throws(
      () => variants.filterRandomly(0.1, { seed: given as number }),
      (error: unknown) => {
        assert.ok(error instanceof Error);
        assert.match(error.message, /`seed` is a whole number from 0 to 9007199254740991/);
        assert.ok(error.message.endsWith(`${written} was given`), error.message);
        return true;
      },
    );
    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

test("filterRandomly with an option it does not know is an Error that names it", () => {
  // A misspelt seed would be a seed that was not given, and the call would
  // keep the sample of 42 with no word of it.
  const variants = many();
  assert.throws(
    () => variants.filterRandomly(0.1, { sed: 7 } as { seed?: number }),
    {
      name: "Error",
      message:
        "popnei: `sed` is not an option of `filterRandomly`, whose options are `seed`",
    },
  );
  assert.deepEqual(variants.steps, []);
  variants.free();
});

test("filterRandomly refuses a second filter of its kind with both keep rates and seeds", () => {
  const variants = many();
  variants.filterRandomly(0.1, { seed: 42 });
  assert.throws(() => variants.filterRandomly(0.5, { seed: 7 }), {
    name: "Error",
    message:
      /filtered by random already, with a keep rate of 0\.1 and a seed of 42, and a second filter of that kind, with a keep rate of 0\.5 and a seed of 7,/,
  });
  assert.deepEqual(variants.steps, [
    { kind: "random", args: { keepRate: 0.1, seed: 42 } },
  ]);
  variants.free();
});

test("filterRandomly after filterFirstN is an Error, as every filter of the variants is", () => {
  const variants = many();
  variants.filterFirstN(10);
  assert.throws(() => variants.filterRandomly(0.1), {
    name: "Error",
    message: /filtered by first_n already, and a filter by random after it/,
  });
  assert.deepEqual(variants.steps, [
    { kind: "first_n", args: { numVars: 10 } },
  ]);
  variants.free();
});
