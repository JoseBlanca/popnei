/**
 * The filter by linkage disequilibrium from TypeScript: which variants
 * `filterByLd` keeps, what it counts, the step it puts on a `Variants` and
 * what it refuses.
 *
 * The item "The filter by linkage disequilibrium" of
 * `docs/specs/filters.md` has the rule and the numbers. The file it is run
 * on is `tests/reference/ld/ld.vcf.gz` of `docs/specs/ld.md`: two
 * chromosomes of 250 biallelic variants each, 1000 bp apart, of 100 diploid
 * individuals with 3 in 100 genotypes missing, 68 of whose 500 variants
 * hold one dosage and are dropped at every threshold. How many variants
 * each setting keeps and the first five it keeps, by position, are the
 * table of "How it is verified" of that item, which was made against the r²
 * matrix that plink2 v2.0.0-a.7.7 wrote for the same file, and they are
 * here as literals. The cargo tests of the core assert the same four rows
 * over blocks of three sizes; what this file adds is that the rule gives
 * them under WebAssembly, where the rows of a block are read one after
 * another and not on the threads of rayon, and that a user reaches it
 * through `filterByLd`.
 *
 * The bytes of the file are read into a `Uint8Array` and given to
 * `openVcf`, which is how a page gives popnei a file: a tab has no
 * filesystem.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Variants } from "popnei";
import { init, openVcf } from "popnei";

import { referenceLd } from "./reference.ts";

await init();

/** The bytes of the dataset, read once for every test that runs on it. */
const LD_VCF = await referenceLd("ld.vcf.gz");

/** How many variants the dataset holds, every one of them given. */
const NUM_VARS = 500;

/**
 * The four rows of the table of "How it is verified": the window in base
 * pairs, the largest r² a kept variant may have against a variant of that
 * window, how many of the 500 variants are kept and the first five of them,
 * by position. Every one of the five is on chr1.
 */
const THE_TABLE: {
  maxDist: number;
  maxAllowedR2: number;
  kept: number;
  firstFive: number[];
}[] = [
  {
    maxDist: 10000,
    maxAllowedR2: 0.1,
    kept: 84,
    firstFive: [1000, 10000, 16000, 22000, 27000],
  },
  {
    maxDist: 10000,
    maxAllowedR2: 0.3,
    kept: 133,
    firstFive: [1000, 5000, 7000, 11000, 15000],
  },
  {
    maxDist: 50000,
    maxAllowedR2: 0.3,
    kept: 85,
    firstFive: [1000, 5000, 7000, 11000, 15000],
  },
  {
    maxDist: 250000,
    maxAllowedR2: 0.3,
    kept: 85,
    firstFive: [1000, 5000, 7000, 11000, 15000],
  },
];

/**
 * The row of the table that the rest of the tests are run at, 133 variants
 * kept of the 500 with a window of 10000 bp and a threshold of 0.3.
 */
const THE_ROW_OF_133 = THE_TABLE[1] as (typeof THE_TABLE)[number];

/**
 * The row of the table whose window, 250000 bp, is the whole of each of the
 * two chromosomes of the file, so that every variant kept on a chromosome
 * is in the window of every variant after it.
 */
const THE_ROW_OF_A_WHOLE_CHROMOSOME = THE_TABLE[3] as (typeof THE_TABLE)[number];

/**
 * The position of the first variant kept on chr2, which is the first
 * variant of that chromosome at every setting of the table: a window ends
 * at the chromosome of the variant it is of, so nothing kept on chr1 is
 * compared with it.
 */
const FIRST_KEPT_OF_CHR2 = 1000;

/**
 * The size of block one test asks for, which makes 72 blocks of the file
 * instead of one: the rule reads the position of a variant and never where
 * a block was cut, so the variants kept are the same.
 */
const NUM_VARS_PER_BLOCK = 7;

/**
 * The numbers that no threshold of a filter takes, and how the message
 * spells each of them: as JavaScript spells a number, `Infinity` and not
 * the `inf` of Rust. A pyNei user writes 0.1 for a threshold on the
 * absolute value of r, which is a `maxAllowedR2` of 0.01 here, and no
 * number outside 0 to 1 is an r² at all.
 */
const THRESHOLDS_REFUSED: [number, string][] = [
  [Number.NaN, "NaN"],
  [-0.1, "-0.1"],
  [1.5, "1.5"],
  [Number.POSITIVE_INFINITY, "Infinity"],
];

/** What is not a number at all, with what the message says was given. */
const NOT_NUMBERS: [unknown, string][] = [
  [undefined, "undefined"],
  [null, "null"],
  ["0.3", "the string `0.3`"],
];

/** The windows that no filter takes, with what the message says was given. */
const WINDOWS_REFUSED: [unknown, string][] = [
  [0, "the number 0"],
  [-10000, "the number -10000"],
  [1500.5, "the number 1500.5"],
  [undefined, "undefined"],
];

/**
 * The largest window a user can ask for, 2^53 - 1 base pairs, and the first
 * one above it.
 *
 * The core takes a window of up to 2^64 - 1, which is what a user of popnei
 * in Python writes. A number of JavaScript counts in twos above 2^53 - 1,
 * so 9007199254740992 written here would not be the window that arrived,
 * and it is refused instead.
 */
const LARGEST_WINDOW = 9007199254740991;
const FIRST_WINDOW_REFUSED = 9007199254740992;

/** The 500 variants of `ld.vcf.gz`, every one of them given. */
function theDataset(): Variants {
  return openVcf(LD_VCF);
}

/**
 * The chromosome and the position of every variant one whole pass over
 * `variants` gives, and the counts of that pass.
 */
function keptBy(
  variants: Variants,
  numVarsPerBlock?: number,
): { chroms: string[]; positions: number[]; passStats: PassStats } {
  const blocks = variants.iterBlocks({
    fields: ["chrom", "pos"],
    ...(numVarsPerBlock === undefined ? {} : { numVarsPerBlock }),
  });
  const chroms: string[] = [];
  const positions: number[] = [];
  for (const block of blocks) {
    if (block.chrom === null || block.pos === null) {
      throw new Error(
        "the pass was asked for the chromosomes and the positions and gave none",
      );
    }
    chroms.push(...block.chrom);
    positions.push(...block.pos);
  }
  return { chroms, positions, passStats: blocks.passStats };
}

for (const { maxDist, maxAllowedR2, kept, firstFive } of THE_TABLE) {
  test(`the filter at an r² of ${maxAllowedR2} within ${maxDist} bp keeps the ${kept} variants of the table`, () => {
    const variants = theDataset();
    variants.filterByLd(maxAllowedR2, maxDist);

    const { chroms, positions, passStats } = keptBy(variants);

    assert.equal(positions.length, kept);
    assert.deepEqual(positions.slice(0, 5), firstFive);
    assert.deepEqual(chroms.slice(0, 5), [
      "chr1",
      "chr1",
      "chr1",
      "chr1",
      "chr1",
    ]);
    assert.equal(positions[chroms.indexOf("chr2")], FIRST_KEPT_OF_CHR2);
    assert.equal(passStats.numVars, kept);
    assert.deepEqual(passStats.filtering, {
      ld: { varsProcessed: NUM_VARS, varsKept: kept },
    });
    variants.free();
  });
}

test("the variants kept do not change with the size of the blocks", () => {
  // Whether a variant is kept turns on the variants kept before it, which
  // the filter carries from one block to the next, so a block that ends in
  // the middle of a window has to change nothing: 72 blocks of 7 variants
  // keep what one block of 500 keeps.
  const variants = theDataset();
  variants.filterByLd(THE_ROW_OF_133.maxAllowedR2, THE_ROW_OF_133.maxDist);

  const { positions, passStats } = keptBy(variants, NUM_VARS_PER_BLOCK);

  assert.equal(positions.length, THE_ROW_OF_133.kept);
  assert.deepEqual(positions.slice(0, 5), THE_ROW_OF_133.firstFive);
  assert.deepEqual(passStats.filtering, {
    ld: { varsProcessed: NUM_VARS, varsKept: THE_ROW_OF_133.kept },
  });
  variants.free();
});

test("filterByLd returns nothing and adds its step with both arguments", () => {
  // The method changes the `Variants` and gives nothing back, as the three
  // threshold filters do, and its step carries the two numbers the user
  // wrote under the names they wrote them in.
  const variants = theDataset();
  assert.deepEqual(variants.steps, []);

  assert.equal(variants.filterByLd(0.3, 10000), undefined);

  assert.deepEqual(variants.steps, [
    { kind: "ld", args: { maxAllowedR2: 0.3, maxDist: 10000 } },
  ]);
  variants.free();
});

test("a maf filter before it is counted apart and is the one given the 500 variants", () => {
  // pyNei filters by the major allele frequency and by linkage
  // disequilibrium in one call with one pair of counts. Here they are two
  // steps: the maf filter reads the source, the filter by linkage
  // disequilibrium reads what it kept, and each has its own counts, in the
  // order of the steps. What the two keep of this file is no number the
  // spec gives, so what is asserted is how the counts fit together.
  const variants = theDataset();
  variants.filterByMaf(0.95);
  variants.filterByLd(0.3, 10000);

  const { positions, passStats } = keptBy(variants);

  assert.deepEqual(Object.keys(passStats.filtering), ["maf", "ld"]);
  const maf = passStats.filtering["maf"];
  const ld = passStats.filtering["ld"];
  if (maf === undefined || ld === undefined) {
    throw new Error("the pass counted neither of the two filters it ran");
  }
  assert.equal(maf.varsProcessed, NUM_VARS);
  assert.equal(ld.varsProcessed, maf.varsKept);
  assert.equal(ld.varsKept, positions.length);
  assert.equal(passStats.numVars, positions.length);
  assert.ok(ld.varsKept < maf.varsKept);
  variants.free();
});

for (const [maxAllowedR2, written] of THRESHOLDS_REFUSED) {
  test(`the filter refuses a maxAllowedR2 of ${written} at the call`, () => {
    // r² is a correlation squared, so no number outside 0 to 1 says which
    // variants a user wants to be rid of. The message names the argument
    // and the value as the user wrote it, and no step is added.
    const variants = theDataset();

    assert.throws(
      () => variants.filterByLd(maxAllowedR2, 10000),
      (error: unknown) =>
        error instanceof Error &&
        error.message.includes(`\`maxAllowedR2\` is ${written},`),
    );

    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

for (const [maxAllowedR2, written] of NOT_NUMBERS) {
  test(`the filter refuses a maxAllowedR2 that is ${written}`, () => {
    // Each of these would reach the core as a number of its own, `null` as
    // a threshold of 0, which keeps one variant of every pair whose
    // dosages say anything about each other at all: what a user gets has
    // to be an error and not a filter they did not write.
    const variants = theDataset();

    assert.throws(
      () => variants.filterByLd(maxAllowedR2 as number, 10000),
      (error: unknown) =>
        error instanceof Error &&
        error.message.includes("maxAllowedR2") &&
        error.message.includes(written),
    );

    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

for (const [maxDist, written] of WINDOWS_REFUSED) {
  test(`the filter refuses a maxDist that is ${written}`, () => {
    // The window is a whole number of base pairs of 1 or more: 0 reaches
    // no variant but the ones at the very position of the variant it is
    // the window of, and a negative number or a fraction would arrive at
    // the core as another number altogether.
    const variants = theDataset();

    assert.throws(
      () => variants.filterByLd(0.3, maxDist as number),
      (error: unknown) =>
        error instanceof Error &&
        error.message.includes("maxDist") &&
        error.message.includes(written),
    );

    assert.deepEqual(variants.steps, []);
    variants.free();
  });
}

test("the largest window a number of JavaScript holds exactly is taken", () => {
  // The core takes the window as a 64 bit whole number of base pairs, and
  // what JavaScript hands it has to be the number the user wrote: at
  // 9007199254740991 it still is. A window longer than a chromosome is the
  // window of the whole chromosome, so this keeps what the row of 250000 bp
  // of the table keeps, and a window that arrived as 0 or as a number cut
  // short would keep another set.
  const variants = theDataset();
  variants.filterByLd(
    THE_ROW_OF_A_WHOLE_CHROMOSOME.maxAllowedR2,
    LARGEST_WINDOW,
  );

  assert.deepEqual(variants.steps, [
    {
      kind: "ld",
      args: {
        maxAllowedR2: THE_ROW_OF_A_WHOLE_CHROMOSOME.maxAllowedR2,
        maxDist: LARGEST_WINDOW,
      },
    },
  ]);
  const { positions } = keptBy(variants);
  assert.equal(positions.length, THE_ROW_OF_A_WHOLE_CHROMOSOME.kept);
  assert.deepEqual(
    positions.slice(0, 5),
    THE_ROW_OF_A_WHOLE_CHROMOSOME.firstFive,
  );
  variants.free();
});

test("a window above the whole numbers a float64 holds is refused", () => {
  // 9007199254740992 is 2^53, where a number of JavaScript starts counting
  // in twos: the core would be given a window the user cannot write and
  // cannot read back, so the package refuses it and says why.
  const variants = theDataset();

  assert.throws(
    () => variants.filterByLd(0.3, FIRST_WINDOW_REFUSED),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("maxDist") &&
      error.message.includes(`${FIRST_WINDOW_REFUSED}`) &&
      error.message.includes("holds exactly"),
  );

  assert.deepEqual(variants.steps, []);
  variants.free();
});

test("a second filter by linkage disequilibrium is refused with the one that is set", () => {
  // Two of them on one `Variants` keep what the stricter of the two keeps
  // alone, so the second says that the steps are not what their user
  // thinks, which running a cell of a notebook twice gives. The message
  // names the kind and both thresholds, and the first step stays as it
  // was.
  const variants = theDataset();
  variants.filterByLd(0.3, 10000);

  assert.throws(
    () => variants.filterByLd(0.1, 50000),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("ld") &&
      error.message.includes("0.3") &&
      error.message.includes("0.1"),
  );

  assert.deepEqual(variants.steps, [
    { kind: "ld", args: { maxAllowedR2: 0.3, maxDist: 10000 } },
  ]);
  variants.free();
});

test("filterByLd throws after the variants were freed", () => {
  // The steps live in the memory of wasm, which `free` gives back, so
  // there is nothing for the filter to be added to.
  const variants = theDataset();
  variants.free();

  assert.throws(() => variants.filterByLd(0.3, 10000), {
    name: "Error",
    message: /freed, so their steps cannot be read or changed/,
  });
});
