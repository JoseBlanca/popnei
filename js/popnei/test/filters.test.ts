/**
 * The three filters from TypeScript: which variants they keep, what they
 * count and what a `Variants` carries once they are put on it.
 *
 * `docs/specs/filters.md` has the three filters, the counts of each and the
 * step that a filter is in a `Variants`. The file they are run on is
 * `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50 diploid
 * individuals, read with every variant given, those that failed their
 * FILTER too, which is what the table of the spec was made on. Its numbers
 * are here as literals: how many variants each filter keeps and the first
 * five it keeps, by position, are what bcftools 1.24 keeps and what pyNei
 * at ef0ca6e keeps, and `tests/reference/filters/` holds every position of
 * each of them. The comparison with pyNei itself is the one of
 * `tests/test_filters.py`, which runs both libraries; node runs neither.
 *
 * This is the first file in which a filter of popnei runs under wasm, where
 * the rows of a block are read one after another and not on the threads of
 * rayon.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Variants } from "popnei";
import { init, openVars, openVcf, writeVars } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The variants of `many.vcf`, read with every variant given. */
const MANY_NUM_VARS = 500;

/**
 * The size of block the test that adds a filter from inside a loop of
 * blocks asks for: 7 variants are 72 blocks of the 500, so the pass goes on
 * for many more blocks after the filter was added to the `Variants`.
 */
const NUM_VARS_PER_BLOCK = 7;

/** Which number of a variant a filter compares, as its counts name it. */
type Kind = "missing_data" | "maf" | "obs_het";

/**
 * The method of each filter and the name of its argument, under the kind
 * the counts of that filter have.
 */
const FILTERS: Record<
  Kind,
  {
    filter: (variants: Variants, threshold: number) => void;
    argument: string;
  }
> = {
  missing_data: {
    filter: (variants, threshold) => variants.filterByMissingData(threshold),
    argument: "maxAllowedMissingRate",
  },
  maf: {
    filter: (variants, threshold) => variants.filterByMaf(threshold),
    argument: "maxAllowedMaf",
  },
  obs_het: {
    filter: (variants, threshold) => variants.filterByObsHet(threshold),
    argument: "maxAllowedObsHet",
  },
};

/**
 * The first row of each filter of the table of "How it is verified" of
 * `docs/specs/filters.md`: the filter with its threshold, how many of the
 * 500 variants of `many.vcf` it keeps and the first five it keeps, by
 * position.
 */
const THE_FIRST_ROWS: [Kind, number, number, number[]][] = [
  ["missing_data", 0, 26, [1259, 2110, 2480, 3072, 3257]],
  ["maf", 0.5, 35, [1074, 1296, 1481, 1962, 2110]],
  ["obs_het", 0.1, 22, [1185, 3516, 3923, 4515, 5921]],
];

/**
 * The chain of "How it is verified" of the counts: the missing data filter
 * at 0.04, the maf filter at 0.8 after it and the observed heterozygosity
 * one at 0.5 after that.
 */
const THE_CHAIN: [Kind, number][] = [
  ["missing_data", 0.04],
  ["maf", 0.8],
  ["obs_het", 0.5],
];

/** What each filter of that chain is given and keeps, in the order of the
 * steps, and the first three variants the chain keeps, by position. */
const THE_COUNTS_OF_THE_CHAIN = {
  missing_data: { varsProcessed: 500, varsKept: 215 },
  maf: { varsProcessed: 215, varsKept: 163 },
  obs_het: { varsProcessed: 163, varsKept: 106 },
};
const VARS_KEPT_BY_THE_CHAIN = 106;
const FIRST_KEPT_BY_THE_CHAIN = [1111, 1407, 1518];

/** What the missing data filter at 0.04 alone is given and keeps. */
const MISSING_DATA_AT_0_04 = {
  missing_data: { varsProcessed: 500, varsKept: 215 },
};
const VARS_KEPT_AT_0_04 = 215;

/**
 * The thresholds that are not a number from 0 to 1, which every filter
 * refuses at the call, and how each one is written in the message. A
 * threshold that is not given arrives as `undefined`, which is the call a
 * user makes when they take the default that pyNei has and popnei does not.
 */
const THRESHOLDS_REFUSED: [number, string][] = [
  [undefined as unknown as number, "undefined"],
  [Number.NaN, "NaN"],
  [-0.1, "-0.1"],
  [1.5, "1.5"],
];

/** The 500 variants of `many.vcf`, the ones that failed their FILTER among
 * them, which is what the table of the spec was made on. */
function many(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/** `variants` with each of `filters`, a kind and a threshold, put on it in
 * that order. */
function filtered(
  variants: Variants,
  filters: readonly [Kind, number][],
): Variants {
  for (const [kind, threshold] of filters) {
    FILTERS[kind].filter(variants, threshold);
  }
  return variants;
}

/** The position of every variant one whole pass over `variants` gives, and
 * the counts of that pass. */
function keptBy(variants: Variants): {
  positions: number[];
  passStats: PassStats;
} {
  const blocks = variants.iterBlocks({ fields: ["pos"] });
  const positions: number[] = [];
  for (const block of blocks) {
    if (block.pos === null) {
      throw new Error("the pass was asked for the positions and gave none");
    }
    positions.push(...block.pos);
  }
  return { positions, passStats: blocks.passStats };
}

for (const [kind, threshold, kept, firstFive] of THE_FIRST_ROWS) {
  test(`the ${kind} filter at ${threshold} keeps the variants of the table`, () => {
    const variants = filtered(many(), [[kind, threshold]]);

    const { positions, passStats } = keptBy(variants);

    assert.equal(positions.length, kept);
    assert.deepEqual(positions.slice(0, 5), firstFive);
    assert.equal(passStats.numVars, kept);
    assert.deepEqual(passStats.filtering, {
      [kind]: { varsProcessed: MANY_NUM_VARS, varsKept: kept },
    });
    variants.free();
  });
}

test("the three filters chained keep 106 variants and count what each was given", () => {
  // The missing data filter at 0.04, the maf filter at 0.8 and the observed
  // heterozygosity one at 0.5, which the chained commands of bcftools keep
  // 215, 163 and 106 variants with. The counts come in the order of the
  // steps, which is the reverse of the chain of readers: the outermost
  // filter, the last step, is the one the core gives first.
  const variants = filtered(many(), THE_CHAIN);

  const { positions, passStats } = keptBy(variants);

  assert.equal(positions.length, VARS_KEPT_BY_THE_CHAIN);
  assert.deepEqual(positions.slice(0, 3), FIRST_KEPT_BY_THE_CHAIN);
  assert.equal(passStats.numVars, VARS_KEPT_BY_THE_CHAIN);
  assert.deepEqual(Object.keys(passStats.filtering), [
    "missing_data",
    "maf",
    "obs_het",
  ]);
  assert.deepEqual(passStats.filtering, THE_COUNTS_OF_THE_CHAIN);
  variants.free();
});

test("the three methods return nothing and add their step in order", () => {
  // Each method changes the `Variants` and gives nothing back, as
  // `array.sort()` does, so `const v2 = v1.filterByMaf(0.95)` gives an
  // `undefined` and an error at the next line instead of two names for one
  // filtered object. The thresholds are under the names of the arguments
  // the user wrote them in.
  const variants = many();
  assert.deepEqual(variants.steps, []);

  assert.equal(variants.filterByMissingData(0.04), undefined);
  assert.equal(variants.filterByMaf(0.8), undefined);
  assert.equal(variants.filterByObsHet(0.5), undefined);

  assert.deepEqual(variants.steps, [
    { kind: "missing_data", args: { maxAllowedMissingRate: 0.04 } },
    { kind: "maf", args: { maxAllowedMaf: 0.8 } },
    { kind: "obs_het", args: { maxAllowedObsHet: 0.5 } },
  ]);
  variants.free();
});

for (const kind of Object.keys(FILTERS) as Kind[]) {
  for (const [threshold, written] of THRESHOLDS_REFUSED) {
    test(`the ${kind} filter refuses a threshold of ${written} at the call`, () => {
      // The number of a variant that the threshold is compared with is one
      // count of the variant divided by another, so no other threshold says
      // anything about which variants a user wants: pyNei takes them, and a
      // 95 written for 0.95 filters nothing there and says nothing. The
      // message names the argument and the value, and no step is added.
      const variants = many();

      assert.throws(
        () => FILTERS[kind].filter(variants, threshold),
        (error: unknown) =>
          error instanceof Error &&
          error.message.includes(FILTERS[kind].argument) &&
          error.message.includes(written),
      );

      assert.deepEqual(variants.steps, []);
      variants.free();
    });
  }
}

test("a second filter of one kind is refused with the threshold that is set", () => {
  // Two threshold filters of one kind keep the variants that the stricter
  // of them keeps alone, so the second says that the user has lost track of
  // what their `Variants` holds, which running the cell of a notebook twice
  // gives. The message names the kind and both thresholds, and a filter of
  // another kind between the two changes nothing.
  const variants = filtered(many(), [["maf", 0.8]]);

  assert.throws(
    () => variants.filterByMaf(0.95),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("maf") &&
      error.message.includes("0.8") &&
      error.message.includes("0.95"),
  );
  assert.deepEqual(variants.steps, [
    { kind: "maf", args: { maxAllowedMaf: 0.8 } },
  ]);

  variants.filterByMissingData(0.04);
  assert.throws(() => variants.filterByMaf(0.5), {
    name: "Error",
    message: /maf/,
  });
  variants.free();
});

test("a filter added after a whole iteration holds in the next one", () => {
  // A step can be added at any time, also after a pass: a user looks at the
  // variants as they are in the file, puts a filter on the same `Variants`
  // and reads them again.
  const variants = many();

  const beforeIt = keptBy(variants);
  assert.equal(beforeIt.positions.length, MANY_NUM_VARS);
  assert.deepEqual(beforeIt.passStats.filtering, {});

  variants.filterByMissingData(0.04);
  const afterIt = keptBy(variants);

  assert.equal(afterIt.positions.length, VARS_KEPT_AT_0_04);
  assert.deepEqual(afterIt.passStats.filtering, MISSING_DATA_AT_0_04);
  variants.free();
});

test("a filter added inside a loop of blocks takes no variant out of that pass", () => {
  // The pass that runs took the steps when it started, so it gives the 500
  // variants of the file and counts no filter, and the pass after it has the
  // filter and its counts.
  const variants = many();
  const blocks = variants.iterBlocks({
    numVarsPerBlock: NUM_VARS_PER_BLOCK,
  });

  let given = 0;
  for (const block of blocks) {
    given += block.numVars;
    if (variants.steps.length === 0) {
      variants.filterByMissingData(0.04);
    }
  }

  assert.equal(given, MANY_NUM_VARS);
  assert.deepEqual(blocks.passStats, {
    numVars: MANY_NUM_VARS,
    filtering: {},
  });

  const after = keptBy(variants);
  assert.equal(after.positions.length, VARS_KEPT_AT_0_04);
  assert.deepEqual(after.passStats.filtering, MISSING_DATA_AT_0_04);
  variants.free();
});

test("writeVars writes the variants the filter kept and counts them", () => {
  // The pass is the core's, which writes one batch for each block, so the
  // count is of what was written; the file read back holds those variants
  // and no step of the `Variants` it was written from.
  const variants = filtered(many(), [["missing_data", 0.04]]);

  const written = writeVars(variants);

  assert.deepEqual(written.passStats, {
    numVars: VARS_KEPT_AT_0_04,
    filtering: MISSING_DATA_AT_0_04,
  });
  const readBack = openVars(written.bytes);
  assert.deepEqual(readBack.steps, []);
  assert.equal(keptBy(readBack).positions.length, VARS_KEPT_AT_0_04);
  readBack.free();
  variants.free();
});

test("the three filters throw after the variants were freed", () => {
  // The steps live in the memory of wasm, which `free` gives back, so a
  // filter has nothing to be added to.
  const variants = many();
  variants.free();

  for (const kind of Object.keys(FILTERS) as Kind[]) {
    assert.throws(() => FILTERS[kind].filter(variants, 0.5), {
      name: "Error",
      message: /freed/,
    });
  }
});
