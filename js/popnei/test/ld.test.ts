/**
 * The matrix of r² of every pair of variants from TypeScript:
 * `calcRogersHuffR2Matrix` and the `R2Matrix` it gives.
 *
 * "How it is verified" of `docs/specs/ld.md` has the numbers. The dataset
 * is `tests/reference/ld/ld.vcf.gz`: two chromosomes of 250 biallelic
 * variants each, 1000 bp apart, of 100 diploid individuals with 3 in 100
 * genotypes missing, and 68 of its 500 variants have no variance. The five
 * r² asserted here are the table of that section, which plink2
 * v2.0.0-a.7.7 gave on 22 September 2026, and they are written as literals.
 * Nothing here computes an expected value with popnei.
 *
 * The bytes of the file are read into a `Uint8Array` and given to
 * `openVcf`, which is how a page gives popnei a file: a tab has no
 * filesystem.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { R2Matrix, Variants } from "popnei";
import { calcRogersHuffR2Matrix, init, openVcf } from "popnei";

import { referenceLd } from "./reference.ts";

await init();

/** How many variants the dataset holds, every one of them given. */
const NUM_VARS = 500;

/** How many of those 500 variants have one dosage and so no r². */
const NUM_VARS_WITH_NO_VARIANCE = 68;

/**
 * The tolerance of an r², relative, which is the one the spec asks the
 * cargo tests for: it is there for a version of plink2 that computes the
 * expression in another order and not for popnei's own rounding. Every one
 * of the five came out of the wasm build under node equal to the literal to
 * the bit on 23 September 2026, and the cargo test of the same file asserts
 * all 93096 pairs against the matrix plink2 wrote.
 */
const TOLERANCE = 1e-12;

/** The bytes of the dataset, read once for every test that runs on it. */
const LD_VCF = await referenceLd("ld.vcf.gz");

/**
 * The five pairs of the table of "How it is verified" of the spec, each
 * with the chromosome and the position of its two variants and the r²
 * plink2 gives: the two closest pairs of the first variant, one 10000 bp
 * away, one at the far end of its chromosome and one on the other
 * chromosome, which has no distance and gets an r² like any other pair.
 */
const THE_PAIRS_OF_THE_TABLE: {
  chromOfA: string;
  posOfA: number;
  chromOfB: string;
  posOfB: number;
  r2: number;
}[] = [
  { chromOfA: "chr1", posOfA: 1000, chromOfB: "chr1", posOfB: 2000, r2: 0.353466669239891 },
  {
    chromOfA: "chr1",
    posOfA: 1000,
    chromOfB: "chr1",
    posOfB: 3000,
    r2: 0.39849991080910563,
  },
  {
    chromOfA: "chr1",
    posOfA: 1000,
    chromOfB: "chr1",
    posOfB: 11000,
    r2: 0.24053784261608957,
  },
  {
    chromOfA: "chr1",
    posOfA: 1000,
    chromOfB: "chr1",
    posOfB: 250000,
    r2: 0.0256751927810228,
  },
  {
    chromOfA: "chr1",
    posOfA: 1000,
    chromOfB: "chr2",
    posOfB: 1000,
    r2: 0.008140034754693937,
  },
];

/** The `Variants` of the dataset, which each test frees when it is done. */
async function theLdDataset(): Promise<Variants> {
  return openVcf(LD_VCF);
}

/**
 * Which row of the matrix is the variant of `chrom` at `pos`, which is how
 * the cargo test of the same table finds them: the rows are the variants in
 * the order the pass gave them, and `chroms` and `poss` say which variant
 * each row is.
 */
function theVariantAt(matrix: R2Matrix, chrom: string, pos: number): number {
  for (let variant = 0; variant < matrix.numVars; variant += 1) {
    if (matrix.chroms[variant] === chrom && matrix.poss[variant] === pos) {
      return variant;
    }
  }
  throw new Error(`no variant of the matrix is ${chrom}:${pos}`);
}

/** The r² of the row `first` and the column `second` of the matrix. */
function theR2OfThePair(
  matrix: R2Matrix,
  first: number,
  second: number,
): number {
  const r2 = matrix.r2[first * matrix.numVars + second];
  if (r2 === undefined) {
    throw new Error(`the matrix holds no cell for the pair ${first}, ${second}`);
  }
  return r2;
}

/** That `found` is `expected` within [`TOLERANCE`], relative. */
function assertTheR2Is(found: number, expected: number, what: string): void {
  assert.ok(
    Math.abs(found - expected) <= TOLERANCE * Math.abs(expected),
    `${what}: the r² is ${found} and not ${expected}`,
  );
}

test("the five r² of the table are the ones plink2 gives", async () => {
  const variants = await theLdDataset();
  try {
    const matrix = calcRogersHuffR2Matrix(variants);
    assert.equal(matrix.numVars, NUM_VARS);
    assert.equal(matrix.r2.length, NUM_VARS * NUM_VARS);
    assert.equal(matrix.chroms.length, NUM_VARS);
    assert.equal(matrix.poss.length, NUM_VARS);
    for (const pair of THE_PAIRS_OF_THE_TABLE) {
      const named =
        `the pair of ${pair.chromOfA}:${pair.posOfA} and ` +
        `${pair.chromOfB}:${pair.posOfB}`;
      const first = theVariantAt(matrix, pair.chromOfA, pair.posOfA);
      const second = theVariantAt(matrix, pair.chromOfB, pair.posOfB);
      assertTheR2Is(theR2OfThePair(matrix, first, second), pair.r2, named);
      // The two cells of a pair hold the same value, which is what says
      // that the row and the column of the matrix were not swapped on the
      // way out of WebAssembly.
      assert.equal(
        theR2OfThePair(matrix, second, first),
        theR2OfThePair(matrix, first, second),
        `${named}: the two cells of the pair`,
      );
    }
  } finally {
    variants.free();
  }
});

test("the chromosomes and the positions are the variants of the file", async () => {
  const variants = await theLdDataset();
  try {
    const matrix = calcRogersHuffR2Matrix(variants);
    // Two chromosomes of 250 variants each, 1000 bp apart, as the spec
    // describes the dataset.
    assert.equal(matrix.chroms[0], "chr1");
    assert.equal(matrix.poss[0], 1000);
    assert.equal(matrix.chroms[249], "chr1");
    assert.equal(matrix.poss[249], 250000);
    assert.equal(matrix.chroms[250], "chr2");
    assert.equal(matrix.poss[250], 1000);
    assert.equal(matrix.chroms[499], "chr2");
    assert.equal(matrix.poss[499], 250000);
    assert.ok(
      matrix.poss instanceof Float64Array,
      "the positions are a Float64Array",
    );
  } finally {
    variants.free();
  }
});

test("a variant of one dosage has NaN in its row, its column and its diagonal", async () => {
  const variants = await theLdDataset();
  try {
    const matrix = calcRogersHuffR2Matrix(variants);
    let withNoVariance = 0;
    for (let variant = 0; variant < matrix.numVars; variant += 1) {
      const diagonal = theR2OfThePair(matrix, variant, variant);
      if (!Number.isNaN(diagonal)) {
        // A variant that has two dosages at least is 1 against itself.
        assert.equal(diagonal, 1, `the variant ${variant} against itself`);
        continue;
      }
      withNoVariance += 1;
      for (let other = 0; other < matrix.numVars; other += 1) {
        assert.ok(
          Number.isNaN(theR2OfThePair(matrix, variant, other)),
          `the row of the variant ${variant} at the column ${other}`,
        );
        assert.ok(
          Number.isNaN(theR2OfThePair(matrix, other, variant)),
          `the column of the variant ${variant} at the row ${other}`,
        );
      }
    }
    assert.equal(withNoVariance, NUM_VARS_WITH_NO_VARIANCE);
  } finally {
    variants.free();
  }
});

test("the counts of the pass are those of the steps of the variants", async () => {
  const variants = await theLdDataset();
  try {
    // A missing data filter of 1 keeps every variant, so the two counts of
    // the filter are the 500 variants of the file and no number here comes
    // from popnei.
    variants.filterByMissingData(1);
    const matrix = calcRogersHuffR2Matrix(variants);
    assert.equal(matrix.numVars, NUM_VARS);
    assert.deepEqual(matrix.passStats, {
      numVars: NUM_VARS,
      filtering: {
        missing_data: { varsProcessed: NUM_VARS, varsKept: NUM_VARS },
      },
    });
  } finally {
    variants.free();
  }
});

test("a pass of more variants than maxNumVars is refused", async () => {
  const variants = await theLdDataset();
  try {
    // The message names the argument of the call a user wrote, which is
    // `maxNumVars` and not the `max_num_vars` the core and the Python
    // package have: what they are told to raise has to be a name of their
    // own language.
    assert.throws(
      () => calcRogersHuffR2Matrix(variants, { maxNumVars: 100 }),
      (error: Error) => {
        assert.ok(
          error.message.includes("`maxNumVars` is 100"),
          `the message is ${error.message}`,
        );
        assert.ok(
          error.message.includes("raise `maxNumVars` or filter the variants"),
          `the message is ${error.message}`,
        );
        assert.ok(
          !error.message.includes("max_num_vars"),
          `the message is ${error.message}`,
        );
        return true;
      },
    );
    // The `Variants` is as it was: the refused pass left nothing behind and
    // the matrix of every variant is still asked for.
    assert.equal(calcRogersHuffR2Matrix(variants).numVars, NUM_VARS);
  } finally {
    variants.free();
  }
});

test("a maxNumVars that is not a whole number of one or more is refused", async () => {
  const variants = await theLdDataset();
  try {
    for (const maxNumVars of [0, -1, 2.5, "500", null]) {
      assert.throws(
        () =>
          calcRogersHuffR2Matrix(variants, {
            maxNumVars: maxNumVars as number,
          }),
        /`maxNumVars` is a whole number of 1 or more/,
        `a maxNumVars of ${String(maxNumVars)}`,
      );
    }
  } finally {
    variants.free();
  }
});

test("a maxNumVars of more variants than the pairs of a browser are counted in is refused", async () => {
  const variants = await theLdDataset();
  try {
    // The matrix holds one value for each pair, the variants squared, and a
    // whole number of the core is 32 bits wide in a browser: 65535 variants
    // are 4294836225 values and 65536 are 4294967296, which is one more
    // than it counts to. A cap of 100000 was what the package invited
    // before this, with a message that said it took 4294967295.
    for (const maxNumVars of [65536, 100000, 4294967295]) {
      assert.throws(
        () => calcRogersHuffR2Matrix(variants, { maxNumVars }),
        (error: Error) => {
          assert.ok(
            error.message.includes(
              "`maxNumVars` is a whole number of 1 or more and at most 65535",
            ),
            `the message is ${error.message}`,
          );
          assert.ok(
            error.message.includes(
              "the matrix of 65536 variants holds more of them than it counts",
            ),
            `the message is ${error.message}`,
          );
          return true;
        },
        `a maxNumVars of ${maxNumVars}`,
      );
    }
    // 65535 is taken: the pass of this file gives its 500 variants, and a
    // cap above the variants of the source changes nothing.
    assert.equal(
      calcRogersHuffR2Matrix(variants, { maxNumVars: 65535 }).numVars,
      NUM_VARS,
    );
  } finally {
    variants.free();
  }
});

test("the names of the chromosomes cannot be written into", async () => {
  const variants = await theLdDataset();
  try {
    const matrix = calcRogersHuffR2Matrix(variants);
    // Python gives a tuple here and read only arrays for the numbers, and
    // the names of a `Distances` and the individuals of a `Variants` are
    // frozen in this package: one package answers the question one way.
    // The module is a module, so it is strict and the assignment throws.
    assert.ok(Object.isFrozen(matrix.chroms));
    assert.throws(() => {
      (matrix.chroms as string[])[0] = "chr99";
    }, TypeError);
    assert.equal(matrix.chroms[0], "chr1");
  } finally {
    variants.free();
  }
});

test("something that is not a Variants is refused", () => {
  assert.throws(
    () => calcRogersHuffR2Matrix({} as unknown as Variants),
    /`variants` is what openVcf or openVars gives/,
  );
});

test("a source with no variant says so", async () => {
  const variants = openVcf(
    new TextEncoder().encode(
      [
        "##fileformat=VCFv4.4",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2",
        "",
      ].join("\n"),
    ),
  );
  try {
    assert.throws(
      () => calcRogersHuffR2Matrix(variants),
      /the source has no variant/,
    );
  } finally {
    variants.free();
  }
});
