/**
 * The counts of one pass, and the steps a `Variants` holds.
 *
 * Every consumer of a `Variants` gives back the counts of the pass it made:
 * how many variants it took, and how many variants each filter of the pass
 * was given and kept. `docs/specs/variant.md` has the `passStats` they come
 * in, `docs/specs/filters.md` the counts of one filter and the step that a
 * filter is, and `docs/specs/io_vars.md` what `writeVars` gives back.
 *
 * No test here puts a filter on a `Variants`, so every `filtering` is empty
 * and every `steps` is an empty array: what a filter counts and what it
 * keeps are asserted by `test/filters.test.ts`. The numbers are those of
 * `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of 50 individuals read
 * with every variant given.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Block, Variants } from "popnei";
import { init, openVars, openVcf, writeVars } from "popnei";

import type { PassCounts } from "../wasm/popnei.js";
import { passStatsOf } from "../dist/variant.js";
import { referenceVcf, vcfOf } from "./reference.ts";

await init();

/** The variants of `many.vcf`, read with every variant given. */
const MANY_NUM_VARS = 500;

/**
 * The size of block the tests that stop half way through a pass ask for,
 * and how many blocks of it they take: three blocks of 7 are 21 of the 500
 * variants.
 */
const NUM_VARS_PER_BLOCK = 7;
const BLOCKS_READ = 3;
const VARS_OF_THE_BLOCKS_READ = 21;

/**
 * How many variants a batch of the vars file the tests read holds, which is
 * not the size of the blocks any of them asks for: the counts are of the
 * blocks a user got and not of what the reader took from the file.
 */
const VARS_NUM_VARS_PER_BLOCK = 100;

/**
 * The first position a number of JavaScript does not hold, 2^53 + 1, which
 * the binding crate refuses once per block: it is how a pass fails at an
 * error of its own and not at one of its reader.
 */
const POSITION_ABOVE_THE_LARGEST = "9007199254740993";

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The same 500 variants as a vars file of batches of 100. */
const MANY_VARS = varsFileOfMany();

/** The bytes of the vars file the tests read `many.vcf` from. */
function varsFileOfMany(): Uint8Array {
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  try {
    return writeVars(variants, {
      numVarsPerBlock: VARS_NUM_VARS_PER_BLOCK,
    }).bytes;
  } finally {
    variants.free();
  }
}

/**
 * The 500 variants of `many.vcf`, from the VCF and from a vars file.
 *
 * Every source gives the same counts: what a pass counts is the variants of
 * the blocks it gave, whatever file they were read from.
 */
const MANY: [string, () => Variants][] = [
  ["a vcf", () => openVcf(MANY_VCF, { onlyPassed: false })],
  ["a vars file", () => openVars(MANY_VARS)],
];

/** The variants of every block of one pass. */
function numVarsOf(blocks: Iterable<Block>): number {
  let given = 0;
  for (const block of blocks) {
    given += block.numVars;
  }
  return given;
}

for (const [source, many] of MANY) {
  test(`a whole iteration of ${source} counts every variant it gave`, () => {
    const variants = many();
    const blocks = variants.iterBlocks();

    assert.equal(numVarsOf(blocks), MANY_NUM_VARS);
    assert.deepEqual(blocks.passStats, {
      numVars: MANY_NUM_VARS,
      filtering: {},
    });
    variants.free();
  });

  test(`the counts of a pass over ${source} that is not over are of the blocks it gave`, () => {
    const variants = many();
    // Three blocks of 7 variants, taken out of a source of 500: the counts
    // are read while the pass runs, which is what a user does in the loop
    // of an `iterBlocks`, and they hold the 21 variants of the three blocks
    // and not the variants the reader has taken out of the file, which for
    // the vars file are its batches of 100.
    const blocks = variants.iterBlocks({
      numVarsPerBlock: NUM_VARS_PER_BLOCK,
    });

    assert.deepEqual(blocks.passStats, { numVars: 0, filtering: {} });
    for (let block = 0; block < BLOCKS_READ; block += 1) {
      assert.equal(blocks.next().value?.numVars, NUM_VARS_PER_BLOCK);
    }

    assert.deepEqual(blocks.passStats, {
      numVars: VARS_OF_THE_BLOCKS_READ,
      filtering: {},
    });
    blocks.return?.();
    variants.free();
  });

  test(`every pass over ${source} counts its own variants`, () => {
    const variants = many();
    // A second iteration gives 500 again and not 1000: every call starts a
    // pass of its own, with counts of its own.
    const first = variants.iterBlocks();
    assert.equal(numVarsOf(first), MANY_NUM_VARS);

    const second = variants.iterBlocks();
    assert.equal(numVarsOf(second), MANY_NUM_VARS);
    assert.deepEqual(first.passStats, {
      numVars: MANY_NUM_VARS,
      filtering: {},
    });
    assert.deepEqual(second.passStats, {
      numVars: MANY_NUM_VARS,
      filtering: {},
    });
    variants.free();
  });

  test(`writeVars of ${source} gives the bytes of the file and the counts of its pass`, () => {
    const variants = many();
    // The pass is the core's, which writes one batch for each block, so the
    // count is what was written and not what a loop of TypeScript saw.
    const written = writeVars(variants, {
      numVarsPerBlock: VARS_NUM_VARS_PER_BLOCK,
    });

    assert.deepEqual(written.passStats, {
      numVars: MANY_NUM_VARS,
      filtering: {},
    });
    const read = openVars(written.bytes);
    assert.equal(numVarsOf(read.iterBlocks()), MANY_NUM_VARS);
    read.free();
    variants.free();
  });

  test(`a Variants of ${source} that was just opened has no step`, () => {
    const variants = many();

    assert.deepEqual(variants.steps, []);
    assert.ok(Array.isArray(variants.steps));
    variants.free();
  });
}

test("the counts of a pass are read after it gave its memory back", () => {
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  // The pass lives in the memory of wasm and is freed when the iteration
  // ends, so what a user reads afterwards are the counts the package took
  // from it before it was freed.
  const blocks = variants.iterBlocks({ numVarsPerBlock: NUM_VARS_PER_BLOCK });
  assert.equal(numVarsOf(blocks), MANY_NUM_VARS);

  assert.deepEqual(blocks.passStats, {
    numVars: MANY_NUM_VARS,
    filtering: {},
  });
  // Read again, they are the same: nothing of the pass is left to add.
  assert.deepEqual(blocks.passStats, {
    numVars: MANY_NUM_VARS,
    filtering: {},
  });
  variants.free();
});

test("a source with no variants counts none", () => {
  // A VCF whose header names three individuals and that has no variant. It
  // is not an error, so the pass is a pass like any other and its count
  // is 0.
  const variants = openVcf(vcfOf([]));

  const blocks = variants.iterBlocks();
  assert.deepEqual([...blocks], []);
  assert.deepEqual(blocks.passStats, { numVars: 0, filtering: {} });

  const written = writeVars(variants);
  assert.deepEqual(written.passStats, { numVars: 0, filtering: {} });
  variants.free();
});

test("the block a pass lost with an error is not among its variants", () => {
  // Three variants read in blocks of one, and a fourth line popnei refuses
  // because it gives one individual four alleles. The block the error
  // happened in never reached the user, so the count is of the three blocks
  // they got and not of the four variants the file holds.
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t40\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0/1/1\t1/1",
    ]),
  );
  const blocks = variants.iterBlocks({ numVarsPerBlock: 1 });

  const given: number[] = [];
  assert.throws(
    () => {
      for (const block of blocks) {
        given.push(block.numVars);
      }
    },
    { name: "Error", message: /ind2/ },
  );

  assert.deepEqual(given, [1, 1, 1]);
  assert.deepEqual(blocks.passStats, { numVars: 3, filtering: {} });
  variants.free();
});

test("the block a pass lost with an error of its own is not among its variants", () => {
  // Three variants read in blocks of one, and a fourth at a position that
  // a number of JavaScript does not hold, 2^53 + 1, which the pass itself
  // refuses after the reader gave it the block. That block never reached
  // the user either, so it is not in the count.
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      `chr1\t${POSITION_ABOVE_THE_LARGEST}\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1`,
    ]),
  );
  const blocks = variants.iterBlocks({ numVarsPerBlock: 1 });

  const given: number[] = [];
  assert.throws(
    () => {
      for (const block of blocks) {
        given.push(block.numVars);
      }
    },
    { name: "Error", message: new RegExp(POSITION_ABOVE_THE_LARGEST) },
  );

  assert.deepEqual(given, [1, 1, 1]);
  assert.deepEqual(blocks.passStats, { numVars: 3, filtering: {} });
  variants.free();
});

test("the counts of a pass are read after a break left the iteration", () => {
  // A pass left with a `break` frees its memory of wasm in the `finally`
  // of its generator, and the counts of the blocks the user got are read
  // after that.
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  const blocks = variants.iterBlocks({ numVarsPerBlock: NUM_VARS_PER_BLOCK });

  let taken = 0;
  for (const block of blocks) {
    taken += block.numVars;
    if (taken === VARS_OF_THE_BLOCKS_READ) {
      break;
    }
  }

  assert.equal(taken, VARS_OF_THE_BLOCKS_READ);
  assert.deepEqual(blocks.passStats, {
    numVars: VARS_OF_THE_BLOCKS_READ,
    filtering: {},
  });
  variants.free();
});

test("the counts of the filters come in the order of the steps", () => {
  // The chain of readers gives the outermost filter first, and a user reads
  // the filters in the order in which they were put on the `Variants`. The
  // numbers are those of `docs/specs/filters.md`, the missing data filter
  // at 0.04 and the maf filter at 0.8 on `many.vcf`: 500 variants given
  // and 215 kept, and then 215 given and 163 kept. The chain has the maf
  // filter first, because it is the outermost.
  const ofTheChain = {
    num_vars: () => 163,
    kinds: () => ["maf", "missing_data"],
    vars_processed: () => Float64Array.from([215, 500]),
    vars_kept: () => Float64Array.from([163, 215]),
    free: () => undefined,
  } as unknown as PassCounts;

  const stats = passStatsOf(ofTheChain);

  assert.deepEqual(Object.keys(stats.filtering), ["missing_data", "maf"]);
  assert.deepEqual(stats.filtering["missing_data"], {
    varsProcessed: 500,
    varsKept: 215,
  });
  assert.deepEqual(stats.filtering["maf"], {
    varsProcessed: 215,
    varsKept: 163,
  });
  assert.equal(stats.numVars, 163);
});

test("the counts of a pass and the steps of a Variants are objects of TypeScript", () => {
  // Nothing a user reads holds memory of wasm: the counts are plain
  // numbers and the steps a plain array, so they answer after the
  // `Variants` was freed and a user keeps them for a report of their own.
  const variants = openVcf(MANY_VCF, { onlyPassed: false });
  const blocks = variants.iterBlocks();
  const given = numVarsOf(blocks);
  const steps = variants.steps;
  variants.free();

  assert.equal(given, MANY_NUM_VARS);
  assert.deepEqual(blocks.passStats, {
    numVars: MANY_NUM_VARS,
    filtering: {},
  });
  assert.deepEqual(steps, []);
});
