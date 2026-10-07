/**
 * The filter of the variants that passed their FILTER from TypeScript,
 * `filterPassed`: which variants it keeps, its counts and its step, and the
 * refusal of a second filter of its kind.
 *
 * `docs/specs/filters.md` has the filter, in "The filter of the variants
 * that passed their FILTER". The file is `many.vcf` of
 * `docs/specs/io_vcf.md`, 500 variants of 50 diploid individuals, read with
 * every variant given, those that failed their FILTER too. The variants it
 * keeps are those of `many.bcftools.tsv`, what bcftools 1.24 printed for
 * `many.vcf`, whose FILTER is `PASS` or a dot: the 475 that `bcftools view
 * -H -f .,PASS many.vcf` gave on 6 October 2026.
 *
 * The refusal of a source with no record of whether its variants passed is
 * tested on `tests/reference/vars/of_1_1.vars`, a vars file of format 1.1
 * that `tests/reference/vars/make_of_1_1.py` writes, which the Python tests
 * read too.
 *
 * `keepsPassed` of a `Variants`, which says beforehand whether the filter can
 * run, is here too, from the paragraph "A `Variants` also has
 * `keeps_passed`" of `docs/specs/variant.md`: true for `many.vcf` and for the
 * vars file written from it, false for `of_1_1.vars` and for a vars file
 * that holds no variant, which popnei writes with the genotypes alone.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { PassStats, Variants } from "popnei";
import { calcPairwiseKosmanDists, init, openVars, openVcf, writeVars } from "popnei";

import { referenceVars, referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The variants of `many.vcf`, read with every variant given. */
const MANY_NUM_VARS = 500;

/** How many of them passed their FILTER, by bcftools 1.24. */
const NUM_THAT_PASSED = 475;

/** The first ten variants that passed, all of chr1, as the spec has them;
 * 1259, between the seventh and the eighth, is the first that failed. */
const THE_FIRST_TEN_THAT_PASSED = [
  1000, 1037, 1074, 1111, 1148, 1185, 1222, 1296, 1333, 1370,
];

/** The chromosome and the position of each variant of `many.bcftools.tsv`
 * whose FILTER is `PASS` or a dot, in the order of the file. */
async function theVariantsThatPassedByBcftools(): Promise<string[]> {
  const text = new TextDecoder().decode(await referenceVcf("many.bcftools.tsv"));
  const passed: string[] = [];
  for (const line of text.split("\n")) {
    if (line === "") {
      continue;
    }
    const [chrom, pos, , , , , filter] = line.split("\t");
    if (filter === "PASS" || filter === ".") {
      passed.push(`${chrom}:${pos}`);
    }
  }
  return passed;
}

const THE_VARIANTS_THAT_PASSED = await theVariantsThatPassedByBcftools();

/** The 500 variants of `many.vcf`, the ones that failed their FILTER among
 * them. */
function many(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/** The chromosome and the position of every variant one whole pass over
 * `variants` gives, in blocks of `numVarsPerBlock`, and the counts of that
 * pass. */
function keptBy(
  variants: Variants,
  numVarsPerBlock?: number,
): { variants: string[]; positions: number[]; passStats: PassStats } {
  const blocks = variants.iterBlocks({ fields: ["chrom", "pos"], numVarsPerBlock });
  const kept: string[] = [];
  const positions: number[] = [];
  for (const block of blocks) {
    if (block.chrom === null || block.pos === null) {
      throw new Error("the pass was asked for the chromosomes and the positions and gave none");
    }
    for (const [row, pos] of block.pos.entries()) {
      kept.push(`${block.chrom[row]}:${pos}`);
      positions.push(pos);
    }
  }
  return { variants: kept, positions, passStats: blocks.passStats };
}

test("the list of the variants that passed by bcftools has the 475 of the spec", () => {
  assert.equal(THE_VARIANTS_THAT_PASSED.length, NUM_THAT_PASSED);
});

test("filterPassed keeps the 475 variants whose FILTER is PASS or a dot, in blocks of 7 and of the default size", () => {
  for (const numVarsPerBlock of [7, undefined]) {
    const variants = many();
    assert.equal(variants.filterPassed(), undefined);

    const { variants: kept, positions, passStats } = keptBy(variants, numVarsPerBlock);

    assert.deepEqual(kept, THE_VARIANTS_THAT_PASSED);
    assert.deepEqual(positions.slice(0, 10), THE_FIRST_TEN_THAT_PASSED);
    assert.deepEqual(passStats, {
      numVars: NUM_THAT_PASSED,
      filtering: {
        passed: { varsProcessed: MANY_NUM_VARS, varsKept: NUM_THAT_PASSED },
      },
      stoppedEarly: false,
    });
    variants.free();
  }
});

test("filterPassed is a step of the kind passed with no argument", () => {
  const variants = many();
  variants.filterPassed();
  assert.deepEqual(variants.steps, [{ kind: "passed", args: {} }]);
  variants.free();
});

test("filterPassed refuses a second filter of its kind and leaves the steps as they were", () => {
  const variants = many();
  variants.filterPassed();
  assert.throws(() => variants.filterPassed(), {
    name: "Error",
    message:
      /filtered by passed already, and a second filter of that kind would keep the same variants as the first/,
  });
  assert.deepEqual(variants.steps, [{ kind: "passed", args: {} }]);
  variants.free();
});

test("filterPassed over a vars file of 1.1, which has no record of whether its variants passed, throws at the first calculation", async () => {
  // The file holds the four variants of `cases.vcf`, the one at 200 that
  // failed its filter q10 among them; the filter refuses it rather than
  // taking that one as passed.
  const variants = openVars(await referenceVars("of_1_1.vars"));
  variants.filterPassed();
  assert.throws(() => calcPairwiseKosmanDists(variants), {
    name: "Error",
    message:
      /the variants hold no record of whether they passed their FILTER, so the filter of the variants that passed cannot run on them: a vars file holds it from format 1.2, written from a VCF/,
  });
  variants.free();
});

test("keepsPassed is true for a VCF, whichever way it was opened, and answers after free", () => {
  const every = many();
  const passedAlone = openVcf(MANY_VCF);
  assert.equal(every.keepsPassed, true);
  assert.equal(passedAlone.keepsPassed, true);
  every.free();
  passedAlone.free();
  assert.equal(every.keepsPassed, true);
  assert.equal(passedAlone.keepsPassed, true);
});

test("keepsPassed is true for a vars file written from a VCF", () => {
  const variants = many();
  const written = openVars(writeVars(variants).bytes);
  assert.equal(written.keepsPassed, true);
  variants.free();
  written.free();
});

test("keepsPassed is false for a vars file of 1.1, and a filter does not change it", async () => {
  const variants = openVars(await referenceVars("of_1_1.vars"));
  assert.equal(variants.keepsPassed, false);
  variants.filterPassed();
  assert.equal(variants.keepsPassed, false);
  variants.free();
  assert.equal(variants.keepsPassed, false);
});

test("keepsPassed is false for a vars file that holds no variant, which popnei writes with the genotypes alone", () => {
  // A MAF of at most 0 keeps the variants with one allele alone, and many.vcf
  // has none: the file is written from a pass of no variant.
  const variants = many();
  variants.filterByMaf(0);
  const written = writeVars(variants);
  assert.equal(written.passStats.numVars, 0);
  const empty = openVars(written.bytes);
  assert.equal(empty.keepsPassed, false);
  variants.free();
  empty.free();
});
