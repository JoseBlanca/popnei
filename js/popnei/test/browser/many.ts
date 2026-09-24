/**
 * The variants of `many.vcf` as the tests under node assert them, and the
 * one function that asserts them over a `Variants` whatever it was opened
 * from.
 *
 * `tests/reference/vcf/many.vcf` is the file of `docs/specs/io_vcf.md`, 500
 * variants of 50 diploid individuals in 117346 bytes, 25 of them failing
 * their FILTER. The cases of this directory open it as the bytes of the
 * file, as a `File` of the page, gzipped as a `File`, and written again as
 * a vars file and opened from a `File` of that. What each of them asserts
 * is that the variants that come out are the same ones, so the assertions
 * are here and each case opens its source and calls them.
 *
 * The numbers are the literals of `test/filter_individuals.test.ts`, which
 * reads its positions and its genotypes from the table of
 * `docs/specs/io_vcf.md`, and the 100 variants of a block are the literal
 * of `test/vcf.test.ts`. Nothing here is computed.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`, and Playwright the ones whose name ends in `.browser.ts`.
 */

import type { Variants } from "../../dist/web.js";

import { assertEqual } from "./assert.ts";

/**
 * The 500 variants of `many.vcf`, the 25 that failed their FILTER among
 * them, and its 50 individuals.
 */
const NUM_VARS = 500;
const NUM_INDIVIDUALS = 50;

/** The first two individuals of the file, in the order it names them. */
const FIRST_INDIVIDUALS = ["ind00", "ind01"];

/**
 * Three individuals of the file in an order that is not the file's, which
 * is what says which column of the genotypes holds which individual.
 */
const THE_THREE = ["ind05", "ind00", "ind49"];

/** How many variants a block of the pass holds. */
const NUM_VARS_PER_BLOCK = 100;

/** The positions of the first five variants of the file. */
const FIRST_POSITIONS = [1000, 1037, 1074, 1111, 1148];

/**
 * The genotypes of `ind05`, `ind00` and `ind49`, in that order, at the
 * variants of the positions 1000 and 1074, which the VCF writes as `1|1`,
 * `1/1`, `1/1` and `0/1`, `2|1`, `1|2`. The three differ at 1074, so that
 * row is the one that would fail if the columns came out in the order of
 * the file.
 */
const GTS_AT_1000 = [1, 1, 1, 1, 1, 1];
const GTS_AT_1074 = [0, 1, 2, 1, 1, 2];

/**
 * That `variants`, opened over whichever source, gives the individuals, the
 * blocks, the positions and the genotypes of `many.vcf`.
 *
 * It reads the individuals of the source, puts the filter of the three
 * individuals on it and makes one pass in blocks of 100 variants, so what
 * it asserts is one whole reading of the file and not its first bytes. The
 * `Variants` is left filtered and unfreed: its caller opened it and frees
 * it.
 *
 * The source is opened with every variant of the file, the 25 that failed
 * their FILTER among them, which is what `onlyPassed` false asks a VCF for
 * and what the test under node these numbers come from reads: the 500
 * below are all of them and not the 475 that passed.
 *
 * @param what What the source is. It goes in front of every message, so
 * that a failure says which of the cases it came from.
 * @throws {Error} When a value is not the one the tests under node assert,
 * and when the pass gives no positions although it was asked for them.
 */
export function assertsTheVariantsOfMany(
  what: string,
  variants: Variants,
): void {
  assertEqual(
    `${what}: the individuals`,
    variants.numIndividuals,
    NUM_INDIVIDUALS,
  );
  assertEqual(
    `${what}: the first two individuals`,
    variants.individuals.slice(0, FIRST_INDIVIDUALS.length),
    FIRST_INDIVIDUALS,
  );

  variants.filterIndividuals(THE_THREE);
  assertEqual(
    `${what}: the individuals the filter kept`,
    variants.individuals,
    THE_THREE,
  );

  const positions: number[] = [];
  const gts: number[] = [];
  let numVarsOfTheFirstBlock = 0;
  for (const block of variants.iterBlocks({
    fields: ["pos"],
    numVarsPerBlock: NUM_VARS_PER_BLOCK,
  })) {
    if (block.pos === null) {
      throw new Error(
        `${what}: the pass was asked for the positions and gave none`,
      );
    }
    if (positions.length === 0) {
      numVarsOfTheFirstBlock = block.numVars;
    }
    positions.push(...block.pos);
    // One allele at a time: a whole block of genotypes spread into the
    // arguments of a call is more of them than a JavaScript engine takes.
    for (const allele of block.gts) {
      gts.push(allele);
    }
  }

  assertEqual(
    `${what}: the variants of the first block`,
    numVarsOfTheFirstBlock,
    NUM_VARS_PER_BLOCK,
  );
  assertEqual(`${what}: the variants of the pass`, positions.length, NUM_VARS);
  assertEqual(
    `${what}: the positions of its first five variants`,
    positions.slice(0, FIRST_POSITIONS.length),
    FIRST_POSITIONS,
  );
  // The genotypes of a block are the variants one after another, and each
  // variant its three individuals one after another, two alleles each.
  const ofAVariant = THE_THREE.length * 2;
  assertEqual(
    `${what}: the genotypes of the three at chr1 1000`,
    gts.slice(0, ofAVariant),
    GTS_AT_1000,
  );
  const at1074 = positions.indexOf(1074);
  if (at1074 < 0) {
    throw new Error(`${what}: the pass gave no variant at chr1 1074`);
  }
  assertEqual(
    `${what}: the genotypes of the three at chr1 1074`,
    gts.slice(at1074 * ofAVariant, (at1074 + 1) * ofAVariant),
    GTS_AT_1074,
  );
}
