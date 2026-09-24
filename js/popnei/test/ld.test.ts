/**
 * The r² between variants from TypeScript: `calcRogersHuffR2Matrix` with
 * the `R2Matrix` it gives, and `calcLdAndDistPerPop` with the bins of
 * distance of each population that it gives.
 *
 * "How it is verified" of `docs/specs/ld.md` has the numbers. The dataset
 * is `tests/reference/ld/ld.vcf.gz`: two chromosomes of 250 biallelic
 * variants each, 1000 bp apart, of 100 diploid individuals with 3 in 100
 * genotypes missing, and 68 of its 500 variants have no variance. The five
 * r² asserted here are the table of that section, which plink2
 * v2.0.0-a.7.7 gave on 22 September 2026, and the bins are its three
 * tables, which `docs/reports/ld-method/bins.py` worked out from the r²
 * plink2 gave on 24 September 2026. All of them are written as literals,
 * and nothing here computes an expected value with popnei.
 *
 * The bytes of the file are read into a `Uint8Array` and given to
 * `openVcf`, which is how a page gives popnei a file: a tab has no
 * filesystem.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type {
  LdAndDistPerPop,
  LdBins,
  LdDecay,
  R2Matrix,
  Variants,
} from "popnei";
import {
  calcLdAndDistPerPop,
  calcRogersHuffR2Matrix,
  init,
  openVcf,
} from "popnei";

import { referenceLd, vcfOf } from "./reference.ts";

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
    // The sentence is the core's, the one every calculation over a pass
    // gives, which `dists.test.ts` asserts whole for the distances.
    assert.throws(
      () => calcRogersHuffR2Matrix(variants),
      /the pass gave no variant and its source holds none/,
    );
  } finally {
    variants.free();
  }
});

/**
 * The ten bins of the first table of "How it is verified" of
 * `docs/specs/ld.md`: the distances from 1 to 250000 base pairs cut into
 * ten of 25000, with the smallest and the largest distance of each, both
 * included.
 */
const THE_BOUNDS_OF_THE_TEN_BINS: [number, number][] = [
  [1, 25000],
  [25001, 50000],
  [50001, 75000],
  [75001, 100000],
  [100001, 125000],
  [125001, 150000],
  [150001, 175000],
  [175001, 200000],
  [200001, 225000],
  [225001, 250000],
];

/**
 * That first table itself, the one population of every one of the 100
 * individuals at a `maxAllowedMaf` of 0.95: for each of the ten bins, how
 * many pairs it holds and the mean of their r².
 *
 * Every value is the one the table prints, which
 * `docs/reports/ld-method/bins.py` worked out from the r² plink2
 * v2.0.0-a.7.7 gives for these individuals and these variants, and which
 * `tests/reference/ld/ld.bins.txt` holds again with the standard
 * deviations. The cargo tests of the core assert all three tables with
 * their standard deviations.
 */
const THE_BINS_OF_EVERY_INDIVIDUAL: [number, number][] = [
  [8744, 0.20767885551844031],
  [7815, 0.07890359176062511],
  [6846, 0.03508441302751711],
  [5962, 0.02056917580626069],
  [5140, 0.015026451851395499],
  [4168, 0.011542104404978385],
  [3308, 0.01145438215819949],
  [2447, 0.012095572873545887],
  [1481, 0.015365254218410632],
  [530, 0.013266303346602112],
];

/**
 * The two populations of the second and the third table of that same part:
 * `pop_a` the individuals `i000` to `i049` of `ld.vcf.gz` and `pop_b`
 * `i050` to `i099`.
 */
const THE_TWO_POPS: Record<string, string[]> = {
  pop_a: Array.from(
    { length: 50 },
    (_unused, which) => `i${String(which).padStart(3, "0")}`,
  ),
  pop_b: Array.from(
    { length: 50 },
    (_unused, which) => `i${String(which + 50).padStart(3, "0")}`,
  ),
};

/**
 * How many pairs each of the ten bins holds for each of those two
 * populations and the mean of their r², from the second and the third
 * table, which are of a `maxAllowedMaf` of 0.8.
 */
const THE_BINS_OF_THE_TWO_POPS: Record<string, [number, number][]> = {
  pop_a: [
    [7394, 0.22226316432228382],
    [6564, 0.09426862122352],
    [5648, 0.04565040582359873],
    [4918, 0.030489796800475328],
    [4304, 0.024750005446480792],
    [3540, 0.02220350204444288],
    [2872, 0.01770136991605654],
    [2137, 0.02022495044871507],
    [1240, 0.01792337843122636],
    [438, 0.020745833685396994],
  ],
  pop_b: [
    [7625, 0.21935192592998345],
    [6779, 0.08778756463719926],
    [5968, 0.0442939962514918],
    [5189, 0.0321335996857634],
    [4473, 0.02676089802517462],
    [3567, 0.02136589829323258],
    [2823, 0.02578147369672746],
    [2086, 0.02294629377272581],
    [1275, 0.022448095913909734],
    [415, 0.016086351215632733],
  ],
};

/**
 * How many of the 500 variants each of those three tables keeps, from the
 * same part: 432 at the `maxAllowedMaf` of 0.95 of the first table, and 396
 * and 402 at the 0.8 of `pop_a` and of `pop_b`, worked out over the
 * individuals of each population alone.
 */
const VARS_AT_THE_MAF_OF_THE_FIRST_TABLE = 432;
const VARS_OF_POP_A = 396;
const VARS_OF_POP_B = 402;

/** The distances and the bins the three tables were run with. */
const OF_THE_TABLES = { minDist: 1, maxDist: 250000, numBins: 10 };

/**
 * The first row of the table of the curve of that same part, for the same
 * one population of every one of the 100 individuals at a `maxAllowedMaf`
 * of 0.95 and over the same pairs as the ten bins above: the fitted ρ per
 * base pair, which is 4Nr, four times the effective size of the population
 * times the recombination per base pair; the fitted curve at a distance of
 * 0; and the distance in base pairs at which it has fallen to half of that.
 *
 * R 4.6.1's `optimize` gave the three on 24 September 2026, over the 46441
 * pairs of that population grouped at the 249 distances they fall at, and
 * `tests/reference/ld/ld.decay.txt` holds them again.
 */
const THE_CURVE_OF_EVERY_INDIVIDUAL = {
  rhoPerBp: 0.00031727347196446889,
  r2AtZero: 0.46198347107438015,
  halfDist: 6810.5712522189806,
};

/**
 * The second and the third row of that table, the curves of `pop_a` and of
 * `pop_b` over the same pairs as their ten bins above, from the same run of
 * R 4.6.1's `optimize`.
 *
 * The two populations have 50 individuals each, so they share the curve's
 * ceiling at a distance of 0 and are told apart by the other two values:
 * their half distances are 7530.10 and 7259.81 base pairs, 3.7 per 100
 * apart. Handing every population the curve of the first one is what that
 * gap is here to catch.
 */
const THE_CURVES_OF_THE_TWO_POPS: Record<
  string,
  { rhoPerBp: number; r2AtZero: number; halfDist: number }
> = {
  pop_a: {
    rhoPerBp: 0.00030068285442483295,
    r2AtZero: 0.46942148760330576,
    halfDist: 7530.1038938711654,
  },
  pop_b: {
    rhoPerBp: 0.00031187790821896646,
    r2AtZero: 0.46942148760330576,
    halfDist: 7259.8060755719744,
  },
};

/**
 * How close the fitted ρ per base pair and the half distance have to be to
 * R's, relative, which is the tolerance "How it is verified" of the spec
 * gives them: both are where a search stopped, and R's two optimisers land
 * 2.1e-9 of themselves apart on this dataset, so it is 480 times their own
 * disagreement.
 *
 * The r² at a distance of 0 is the curve's ceiling, which the individuals of
 * the population fix on their own with no search, and the spec compares it
 * within the 1e-12 of [`TOLERANCE`].
 */
const TOLERANCE_OF_THE_FIT = 1e-6;

/** That the mean r² `found` is `expected` within [`TOLERANCE`], relative. */
function assertTheMeanIs(found: number, expected: number, what: string): void {
  assert.ok(
    Math.abs(found - expected) <= TOLERANCE * Math.abs(expected),
    `${what}: the mean r² is ${found} and not ${expected}`,
  );
}

/** The bins of the population `pop`, which the pass has to have given. */
function theBinsOf(ofThePass: LdAndDistPerPop, pop: string): LdBins {
  const bins = ofThePass.perPop[pop];
  if (bins === undefined) {
    throw new Error(`the pass gave no bins for the population ${pop}`);
  }
  return bins;
}

/** The curve of the population `pop`, which the pass has to have given. */
function theCurveOf(ofThePass: LdAndDistPerPop, pop: string): LdDecay {
  const curve = ofThePass.decayPerPop[pop];
  if (curve === undefined) {
    throw new Error(`the pass gave no curve for the population ${pop}`);
  }
  return curve;
}

/** That `found` is `expected` within `tolerance`, relative. */
function assertTheFitIs(
  found: number,
  expected: number,
  tolerance: number,
  what: string,
): void {
  assert.ok(
    Math.abs(found - expected) <= tolerance * Math.abs(expected),
    `${what}: it is ${found} and not ${expected}`,
  );
}

test("the ten bins of one population are the ones plink2 gives", async () => {
  const variants = await theLdDataset();
  try {
    const ofThePass = calcLdAndDistPerPop(variants, {
      ...OF_THE_TABLES,
      maxAllowedMaf: 0.95,
    });

    // With no `pops` there is one population of every individual, named as
    // pyNei names it.
    assert.deepEqual(Object.keys(ofThePass.perPop), ["pop"]);
    assert.deepEqual(ofThePass.numVarsPerPop, {
      pop: VARS_AT_THE_MAF_OF_THE_FIRST_TABLE,
    });
    // The pass counted every variant of the file: the major allele
    // frequency takes variants out of a population and not out of the pass.
    assert.equal(ofThePass.passStats.numVars, NUM_VARS);

    const bins = theBinsOf(ofThePass, "pop");
    for (const values of [
      bins.smallestDist,
      bins.largestDist,
      bins.numPairs,
      bins.meanR2,
      bins.sdR2,
    ]) {
      assert.ok(values instanceof Float64Array, "the bins are Float64Arrays");
      assert.equal(values.length, OF_THE_TABLES.numBins);
    }
    assert.deepEqual(
      [...bins.smallestDist],
      THE_BOUNDS_OF_THE_TEN_BINS.map(([smallest]) => smallest),
    );
    assert.deepEqual(
      [...bins.largestDist],
      THE_BOUNDS_OF_THE_TEN_BINS.map(([, largest]) => largest),
    );
    assert.deepEqual(
      [...bins.numPairs],
      THE_BINS_OF_EVERY_INDIVIDUAL.map(([numPairs]) => numPairs),
    );
    for (const [bin, [, meanR2]] of THE_BINS_OF_EVERY_INDIVIDUAL.entries()) {
      assertTheMeanIs(
        bins.meanR2[bin] as number,
        meanR2,
        `the bin ${bin} of the ten`,
      );
      // A bin of this table holds hundreds of pairs at least, so its
      // standard deviation is a number and not the NaN of a bin with none.
      // The cargo tests of the core assert its value.
      assert.ok(
        Number.isFinite(bins.sdR2[bin]),
        `the bin ${bin} of the ten has a standard deviation`,
      );
    }
  } finally {
    variants.free();
  }
});

test("two populations of one pass count their own pairs and their own variants", async () => {
  const variants = await theLdDataset();
  try {
    const ofThePass = calcLdAndDistPerPop(variants, {
      ...OF_THE_TABLES,
      pops: THE_TWO_POPS,
      maxAllowedMaf: 0.8,
    });

    // The populations come back in the order of the keys of `pops`.
    assert.deepEqual(Object.keys(ofThePass.perPop), ["pop_a", "pop_b"]);
    // At 0.95 both populations keep the same 432 variants, so 0.8 is the
    // threshold that fails when the major allele frequency is worked out
    // over all the individuals instead of over those of the population.
    assert.deepEqual(ofThePass.numVarsPerPop, {
      pop_a: VARS_OF_POP_A,
      pop_b: VARS_OF_POP_B,
    });
    for (const [pop, rows] of Object.entries(THE_BINS_OF_THE_TWO_POPS)) {
      const bins = theBinsOf(ofThePass, pop);
      assert.deepEqual(
        [...bins.numPairs],
        rows.map(([numPairs]) => numPairs),
        `the pairs of ${pop}`,
      );
      for (const [bin, [, meanR2]] of rows.entries()) {
        assertTheMeanIs(
          bins.meanR2[bin] as number,
          meanR2,
          `the bin ${bin} of ${pop}`,
        );
      }
    }
    // Each population gets the curve fitted to its own pairs. The three
    // values of a curve cross from the core in one array of a value for
    // each population, and reading the first population's value for every
    // one of them would leave `pop_b` with `pop_a`'s half distance, 3.7 per
    // 100 away from its own.
    for (const [pop, fitted] of Object.entries(THE_CURVES_OF_THE_TWO_POPS)) {
      const curve = theCurveOf(ofThePass, pop);
      assertTheFitIs(
        curve.rhoPerBp,
        fitted.rhoPerBp,
        TOLERANCE_OF_THE_FIT,
        `the ρ per base pair of ${pop}`,
      );
      assertTheFitIs(
        curve.r2AtZero,
        fitted.r2AtZero,
        TOLERANCE,
        `the r² at a distance of 0 of ${pop}`,
      );
      assertTheFitIs(
        curve.halfDist,
        fitted.halfDist,
        TOLERANCE_OF_THE_FIT,
        `the half distance of ${pop}`,
      );
    }
  } finally {
    variants.free();
  }
});

test("the curve of one population is the one R fits to its pairs", async () => {
  const variants = await theLdDataset();
  try {
    const ofThePass = calcLdAndDistPerPop(variants, {
      ...OF_THE_TABLES,
      maxAllowedMaf: 0.95,
    });

    assert.deepEqual(Object.keys(ofThePass.decayPerPop), ["pop"]);
    const curve = theCurveOf(ofThePass, "pop");
    // The curve of a population is three numbers and not three arrays: one
    // curve is fitted to every pair it counted, at the distance of each
    // pair, so the bins do not cut it and `numBins` does not move it.
    for (const value of [curve.rhoPerBp, curve.r2AtZero, curve.halfDist]) {
      assert.equal(typeof value, "number");
    }
    assertTheFitIs(
      curve.rhoPerBp,
      THE_CURVE_OF_EVERY_INDIVIDUAL.rhoPerBp,
      TOLERANCE_OF_THE_FIT,
      "the ρ per base pair",
    );
    assertTheFitIs(
      curve.r2AtZero,
      THE_CURVE_OF_EVERY_INDIVIDUAL.r2AtZero,
      TOLERANCE,
      "the r² at a distance of 0",
    );
    assertTheFitIs(
      curve.halfDist,
      THE_CURVE_OF_EVERY_INDIVIDUAL.halfDist,
      TOLERANCE_OF_THE_FIT,
      "the half distance",
    );
    // Fitting the mean of each of these ten bins at the middle of the bin
    // instead gives a half distance of 7886.60 bp, 15.8 per 100 above the
    // 6810.57 the pairs give, so a fit that read the bins would fail this.
    assert.ok(
      curve.halfDist < 7000,
      `the half distance is ${curve.halfDist} and the pairs give 6810.57 bp`,
    );
  } finally {
    variants.free();
  }
});

test("a population whose pairs fall at one distance and one with no pair get no curve", async () => {
  // Two populations of the three individuals of the file, one for each way
  // a pass reaches the three NaN of "The cases" of `docs/specs/ld.md`.
  // `of_the_three` keeps every variant, and the `maxDist` of 15 leaves it
  // the two pairs 10 base pairs apart and drops the one 20 apart, so every
  // pair it counts is at one distance, which says nothing about a fall-off.
  // `of_one_individual` keeps the two variants `ind3` is heterozygous at
  // and counts no pair at all: one individual has one dosage at every
  // variant, so no variant of it has variance and its one pair has no r².
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1",
      "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/1",
    ]),
  );
  try {
    const ofThePass = calcLdAndDistPerPop(variants, {
      pops: {
        of_the_three: ["ind1", "ind2", "ind3"],
        of_one_individual: ["ind3"],
      },
      minDist: 1,
      maxDist: 15,
      numBins: 3,
    });

    // The two come back in the order they were given and not in the order
    // of their names, as the bins and the variants of each do.
    assert.deepEqual(Object.keys(ofThePass.decayPerPop), [
      "of_the_three",
      "of_one_individual",
    ]);
    assert.deepEqual(ofThePass.numVarsPerPop, {
      of_the_three: 3,
      of_one_individual: 2,
    });
    // The bins say which of the two cases each population is: the pairs of
    // the first are the two 10 base pairs apart, both in the second bin, of
    // 6 to 10, and the second population has no pair in any bin.
    assert.deepEqual(
      [...theBinsOf(ofThePass, "of_the_three").numPairs],
      [0, 2, 0],
    );
    assert.deepEqual(
      [...theBinsOf(ofThePass, "of_one_individual").numPairs],
      [0, 0, 0],
    );
    for (const pop of Object.keys(ofThePass.decayPerPop)) {
      const curve = theCurveOf(ofThePass, pop);
      assert.ok(Number.isNaN(curve.rhoPerBp), `the ρ per base pair of ${pop}`);
      assert.ok(Number.isNaN(curve.r2AtZero), `the r² at 0 of ${pop}`);
      assert.ok(Number.isNaN(curve.halfDist), `the half distance of ${pop}`);
    }
  } finally {
    variants.free();
  }
});

test("a bin that no pair reaches holds no pair and has no mean and no standard deviation", async () => {
  // Three variants 10 base pairs apart, which no bin from 1000 to 5000
  // reaches: "The cases" of `docs/specs/ld.md` says that every bin empty is
  // no error, and the three variants still passed the major allele
  // frequency, so `numVarsPerPop` counts them.
  const variants = openVcf(
    vcfOf([
      "chr1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t20\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
      "chr1\t30\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/1\t0/0",
    ]),
  );
  try {
    const ofThePass = calcLdAndDistPerPop(variants, {
      minDist: 1000,
      maxDist: 5000,
      numBins: 4,
    });

    assert.deepEqual(ofThePass.numVarsPerPop, { pop: 3 });
    const bins = theBinsOf(ofThePass, "pop");
    assert.deepEqual([...bins.smallestDist], [1000, 2001, 3001, 4001]);
    assert.deepEqual([...bins.largestDist], [2000, 3000, 4000, 5000]);
    assert.deepEqual([...bins.numPairs], [0, 0, 0, 0]);
    assert.ok(
      [...bins.meanR2].every((mean) => Number.isNaN(mean)),
      "every mean of a bin with no pair is NaN",
    );
    assert.ok(
      [...bins.sdR2].every((sd) => Number.isNaN(sd)),
      "every standard deviation of a bin with no pair is NaN",
    );
  } finally {
    variants.free();
  }
});

test("a distance, a number of bins and a frequency that are no number of their own are refused", async () => {
  const variants = await theLdDataset();
  try {
    // A negative distance cannot reach the core, whose distances are
    // unsigned, and neither can one with a fraction: the code wasm-bindgen
    // generates would hand on another number and say nothing. The package
    // refuses both, under the name the user wrote them in.
    for (const minDist of [-1, 2.5, "500", null]) {
      assert.throws(
        () => calcLdAndDistPerPop(variants, { minDist: minDist as number }),
        /`minDist` is a whole number of base pairs of 0 or more/,
        `a minDist of ${String(minDist)}`,
      );
    }
    for (const maxDist of [-250, 2.5, null]) {
      assert.throws(
        () => calcLdAndDistPerPop(variants, { maxDist: maxDist as number }),
        /`maxDist` is a whole number of base pairs of 0 or more/,
        `a maxDist of ${String(maxDist)}`,
      );
    }
    for (const numBins of [-3, 2.5, "ten", null]) {
      assert.throws(
        () => calcLdAndDistPerPop(variants, { numBins: numBins as number }),
        /`numBins` is a whole number of 0 or more/,
        `a numBins of ${String(numBins)}`,
      );
    }
    assert.throws(
      () =>
        calcLdAndDistPerPop(variants, {
          maxAllowedMaf: "a half" as unknown as number,
        }),
      /`maxAllowedMaf` is a number/,
    );
  } finally {
    variants.free();
  }
});

test("an empty range, no bin and a frequency out of 0 to 1 are refused under the names of TypeScript", async () => {
  const variants = await theLdDataset();
  try {
    // The three are the core's rules, and the message a user reads names
    // the argument of the call they wrote and not the `min_dist` of Rust
    // and of Python.
    assert.throws(
      () => calcLdAndDistPerPop(variants, { minDist: 5000, maxDist: 4000 }),
      (error: Error) => {
        assert.ok(
          error.message.includes("`minDist` is 5000 and `maxDist` is 4000"),
          `the message is ${error.message}`,
        );
        assert.ok(
          !error.message.includes("min_dist"),
          `the message is ${error.message}`,
        );
        return true;
      },
    );
    assert.throws(
      () => calcLdAndDistPerPop(variants, { numBins: 0 }),
      (error: Error) => {
        assert.ok(
          error.message.includes("`numBins` is 0"),
          `the message is ${error.message}`,
        );
        assert.ok(
          !error.message.includes("num_bins"),
          `the message is ${error.message}`,
        );
        // The message of `numBins` names the two distances the bins are
        // cut across, and they are the names of this package as well.
        assert.ok(
          !error.message.includes("min_dist"),
          `the message is ${error.message}`,
        );
        assert.ok(
          !error.message.includes("max_dist"),
          `the message is ${error.message}`,
        );
        return true;
      },
    );
    assert.throws(
      () => calcLdAndDistPerPop(variants, { maxAllowedMaf: 1.5 }),
      (error: Error) => {
        assert.ok(
          error.message.includes("`maxAllowedMaf` is 1.5"),
          `the message is ${error.message}`,
        );
        assert.ok(
          !error.message.includes("max_allowed_maf"),
          `the message is ${error.message}`,
        );
        return true;
      },
    );
    // A population that names an individual the pass does not give, which
    // only the pass knows: the refusal is the core's and names it.
    assert.throws(
      () => calcLdAndDistPerPop(variants, { pops: { pop1: ["i000", "i999"] } }),
      /i999/,
    );
    // The `Variants` is as it was: no refused call left anything behind.
    assert.equal(
      calcLdAndDistPerPop(variants, OF_THE_TABLES).passStats.numVars,
      NUM_VARS,
    );
  } finally {
    variants.free();
  }
});

test("the fall-off refuses what is not a Variants as the matrix does", () => {
  assert.throws(
    () => calcLdAndDistPerPop({} as unknown as Variants),
    /`variants` is what openVcf or openVars gives/,
  );
});

test("a call that names no argument counts what the defaults of the core say", async () => {
  const variants = await theLdDataset();
  try {
    const ofThePass = calcLdAndDistPerPop(variants);

    const bins = theBinsOf(ofThePass, "pop");
    // The four defaults are the core's, and what they are is written in
    // "Its Python function" of `docs/specs/ld.md`: the distances from 1 to
    // 1000000 base pairs, cut into 50 bins of 20000, and a `maxAllowedMaf`
    // of 0.95, which 432 of the 500 variants of this file pass.
    assert.equal(bins.numPairs.length, 50);
    assert.equal(bins.smallestDist[0], 1);
    assert.equal(bins.largestDist[0], 20000);
    assert.equal(bins.largestDist[49], 1000000);
    assert.deepEqual(ofThePass.numVarsPerPop, {
      pop: VARS_AT_THE_MAF_OF_THE_FIRST_TABLE,
    });
  } finally {
    variants.free();
  }
});
