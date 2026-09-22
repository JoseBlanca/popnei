/**
 * The filter of individuals from TypeScript: whose genotypes come out of a
 * pass, in which order, and what a `Variants` carries once it is put on.
 *
 * `docs/specs/filters.md` has it under "The filter of individuals". The
 * file it runs on is `many.vcf` of `docs/specs/io_vcf.md`, 500 variants of
 * 50 diploid individuals, read with every variant given, those that failed
 * their FILTER too, which is what pyNei reads and what the numbers of the
 * spec were made on. Those numbers are here as literals: `bcftools view -s
 * ind05,ind00,ind49` keeps the three individuals in that order, and the
 * missing data filter at 0 after it keeps 423 of the 500 variants, where
 * the same filter over the 50 individuals keeps 26. The comparison with
 * pyNei itself is the one of `tests/test_filter_individuals.py`, which runs
 * both libraries; node runs neither.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import type { Variants } from "popnei";
import { init, openVcf } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/** The bytes of `many.vcf`, read once for the whole file. */
const MANY_VCF = await referenceVcf("many.vcf");

/** The 500 variants of `many.vcf`, and its 50 individuals. */
const MANY_NUM_VARS = 500;
const MANY_NUM_INDIVIDUALS = 50;

/**
 * The three individuals of "How it is verified" of the filter, in the order
 * a user names them, which is not the order `many.vcf` has them in.
 */
const THE_THREE = ["ind05", "ind00", "ind49"];

/**
 * What the missing data filter at 0 keeps after the three individuals are
 * taken: 423 variants, at these five positions first. The same filter over
 * the 50 individuals of the source keeps 26.
 */
const KEPT_OF_THE_THREE = 423;
const FIRST_POSITIONS_OF_THE_THREE = [1000, 1037, 1074, 1111, 1148];
const KEPT_OF_THE_FIFTY = 26;

/**
 * The genotypes of `ind05`, `ind00` and `ind49`, in that order, at the
 * variants of the positions 1000 and 1074, which are `1|1`, `1/1`, `1/1`
 * and `0/1`, `2|1`, `1|2` in the VCF. The three differ at 1074, so that row
 * says which column holds which individual: the order of the source would
 * put `ind00` first.
 */
const GTS_AT_1000 = [1, 1, 1, 1, 1, 1];
const GTS_AT_1074 = [0, 1, 2, 1, 1, 2];

/** The 500 variants of `many.vcf`, the ones that failed their FILTER among
 * them, which is what the numbers of the spec were made on. */
function many(): Variants {
  return openVcf(MANY_VCF, { onlyPassed: false });
}

/**
 * The genotypes and the positions of every block of one whole pass over
 * `variants`, joined, with how many individuals its blocks hold.
 */
function joined(variants: Variants): {
  gts: number[];
  positions: number[];
  numIndividuals: number;
} {
  const gts: number[] = [];
  const positions: number[] = [];
  let numIndividuals = 0;
  for (const block of variants.iterBlocks({ fields: ["pos"] })) {
    if (block.pos === null) {
      throw new Error("the pass was asked for the positions and gave none");
    }
    // One allele at a time: a whole block of genotypes spread into the
    // arguments of a call is more of them than a JavaScript engine takes.
    for (const allele of block.gts) {
      gts.push(allele);
    }
    positions.push(...block.pos);
    numIndividuals = block.numIndividuals;
  }
  return { gts, positions, numIndividuals };
}

test("the three individuals come out at every variant in the order they were named", () => {
  // The filter takes columns of the genotypes away and no row, so the 500
  // variants of the file all come out, each with the three genotypes.
  const variants = many();

  variants.filterIndividuals(THE_THREE);
  const { gts, positions, numIndividuals } = joined(variants);

  assert.deepEqual(variants.individuals, THE_THREE);
  assert.equal(numIndividuals, 3);
  assert.equal(positions.length, MANY_NUM_VARS);
  assert.equal(gts.length, MANY_NUM_VARS * 3 * 2);
  assert.deepEqual(gts.slice(0, 6), GTS_AT_1000);
  const at1074 = positions.indexOf(1074);
  assert.deepEqual(gts.slice(at1074 * 6, at1074 * 6 + 6), GTS_AT_1074);
  variants.free();
});

test("the missing data filter after it keeps the 423 variants and counts them", () => {
  // Of the 500 variants, 423 have every genotype of the three called, where
  // only 26 have every genotype of the 50 called: the filter after the step
  // divides by the kept individuals alone. The filter of individuals takes
  // no variant away and has no counts, so the counts of the pass hold the
  // missing data filter alone.
  const variants = many();
  variants.filterIndividuals(THE_THREE);
  variants.filterByMissingData(0);

  const blocks = variants.iterBlocks({ fields: ["pos"] });
  const positions: number[] = [];
  for (const block of blocks) {
    positions.push(...(block.pos ?? []));
  }

  assert.equal(positions.length, KEPT_OF_THE_THREE);
  assert.deepEqual(positions.slice(0, 5), FIRST_POSITIONS_OF_THE_THREE);
  assert.deepEqual(blocks.passStats, {
    numVars: KEPT_OF_THE_THREE,
    filtering: {
      missing_data: {
        varsProcessed: MANY_NUM_VARS,
        varsKept: KEPT_OF_THE_THREE,
      },
    },
  });
  variants.free();
});

test("a filter before it counts over every individual of the source", () => {
  // A step sees what the steps before it gave, so this filter divides by
  // the 50 individuals of the source and keeps the 26 variants with every
  // genotype called, where the same filter after the step keeps 423. The
  // blocks that come out of it hold the three individuals all the same.
  const variants = many();
  variants.filterByMissingData(0);
  variants.filterIndividuals(THE_THREE);

  const { positions, numIndividuals } = joined(variants);

  assert.equal(positions.length, KEPT_OF_THE_FIFTY);
  assert.equal(numIndividuals, 3);
  variants.free();
});

test("the method returns nothing and adds its step with the names it was given", () => {
  // The method changes the `Variants` and gives nothing back, as the three
  // threshold filters do, so `const v2 = v1.filterIndividuals(names)` gives
  // an `undefined` and an error at the next line instead of two names for
  // one filtered object. The names are under the argument the user wrote
  // them in, in their order.
  const variants = many();
  assert.deepEqual(variants.steps, []);

  assert.equal(variants.filterIndividuals(THE_THREE), undefined);

  assert.deepEqual(variants.steps, [
    { kind: "individuals", args: { individuals: THE_THREE } },
  ]);
  assert.deepEqual(variants.iterBlocks().passStats.filtering, {});
  variants.free();
});

test("a threshold filter beside it keeps its threshold in the steps", () => {
  // The arguments of the two kinds of step cross together, a threshold as
  // one number and the individuals as their names, so a `Variants` that
  // holds both is what says that neither is read as the other.
  const variants = many();

  variants.filterByMaf(0.95);
  variants.filterIndividuals(THE_THREE);
  variants.filterByObsHet(0.5);

  assert.deepEqual(variants.steps, [
    { kind: "maf", args: { maxAllowedMaf: 0.95 } },
    { kind: "individuals", args: { individuals: THE_THREE } },
    { kind: "obs_het", args: { maxAllowedObsHet: 0.5 } },
  ]);
  variants.free();
});

test("the individuals are the kept ones before a pass and after one", () => {
  // They are what the next pass gives and not what the header holds, and
  // nothing of a pass changes them.
  const variants = many();
  assert.equal(variants.numIndividuals, MANY_NUM_INDIVIDUALS);
  assert.deepEqual(variants.individuals.slice(0, 2), ["ind00", "ind01"]);

  variants.filterIndividuals(THE_THREE);

  assert.deepEqual(variants.individuals, THE_THREE);
  assert.equal(variants.numIndividuals, 3);

  for (const block of variants.iterBlocks()) {
    assert.equal(block.numIndividuals, 3);
  }

  assert.deepEqual(variants.individuals, THE_THREE);
  assert.equal(variants.numIndividuals, 3);
  variants.free();
});

test("the individuals answer after the variants were freed", () => {
  // The names are in JavaScript and not in the memory of wasm, which is
  // what lets a user read which individuals a result was calculated over
  // after they gave the file back.
  const variants = many();
  variants.filterIndividuals(THE_THREE);

  variants.free();

  assert.deepEqual(variants.individuals, THE_THREE);
  assert.equal(variants.numIndividuals, 3);
});

test("two handles over one file each keep their own individuals", () => {
  // A second filter of individuals on one `Variants` is refused, so a user
  // who wants two sets of individuals opens the source twice, which reads
  // the header and nothing else. Each handle keeps its own, and the
  // genotypes of each are the columns of the whole file at those places.
  const first = many();
  const second = many();
  first.filterIndividuals(["ind00", "ind01"]);
  second.filterIndividuals(["ind02", "ind03", "ind04"]);

  const unfiltered = many();
  const ofTheFirst = joined(first).gts;
  const ofTheSecond = joined(second).gts;
  const whole = joined(unfiltered).gts;
  unfiltered.free();

  assert.deepEqual(first.individuals, ["ind00", "ind01"]);
  assert.deepEqual(second.individuals, ["ind02", "ind03", "ind04"]);
  for (let variant = 0; variant < MANY_NUM_VARS; variant += 1) {
    const row = variant * MANY_NUM_INDIVIDUALS * 2;
    assert.deepEqual(
      ofTheFirst.slice(variant * 4, variant * 4 + 4),
      whole.slice(row, row + 4),
    );
    assert.deepEqual(
      ofTheSecond.slice(variant * 6, variant * 6 + 6),
      whole.slice(row + 4, row + 10),
    );
  }
  first.free();
  second.free();
});

test("a name that is not an individual of the source is refused at the call", () => {
  // pyNei's `filter_samples` drops it in silence and gives a `Variants` of
  // one individual, so a typed name is a result over the wrong individuals
  // there. The message names what the user wrote, the step is not added and
  // the individuals are those of the source.
  const variants = many();

  assert.throws(() => variants.filterIndividuals(["ind05", "nope"]), {
    name: "Error",
    message: /nope/,
  });

  assert.deepEqual(variants.steps, []);
  assert.equal(variants.numIndividuals, MANY_NUM_INDIVIDUALS);
  variants.free();
});

test("a name that is there twice is refused at the call", () => {
  // Two columns of the genotypes of one individual are one individual for
  // everything that reads them, and every count over them would hold it
  // twice. pyNei keeps the individual once.
  const variants = many();

  assert.throws(() => variants.filterIndividuals(["ind05", "ind00", "ind05"]), {
    name: "Error",
    message: /ind05/,
  });

  assert.deepEqual(variants.steps, []);
  variants.free();
});

test("a filter of no individual is refused at the call", () => {
  // It would leave variants of nobody, and every source of popnei holds one
  // individual at least.
  const variants = many();

  assert.throws(() => variants.filterIndividuals([]), { name: "Error" });

  assert.deepEqual(variants.steps, []);
  assert.equal(variants.numIndividuals, MANY_NUM_INDIVIDUALS);
  variants.free();
});

test("a second filter of individuals is refused with the kind", () => {
  // Two lists keep the individuals that are in both, which is one list, so
  // the second says that the user has lost track of the individuals their
  // variants carry, which running the cell of a notebook twice gives. pyNei
  // takes it. After the refusal the steps are as they were, and a threshold
  // filter between the two changes nothing.
  const variants = many();
  variants.filterIndividuals(THE_THREE);

  assert.throws(() => variants.filterIndividuals(["ind01", "ind02"]), {
    name: "Error",
    message: /individuals/,
  });
  assert.deepEqual(variants.steps, [
    { kind: "individuals", args: { individuals: THE_THREE } },
  ]);
  assert.deepEqual(variants.individuals, THE_THREE);

  variants.filterByMaf(0.95);
  assert.throws(() => variants.filterIndividuals(["ind01"]), {
    name: "Error",
    message: /individuals/,
  });
  variants.free();
});

test("one name written as a string is refused before the call", () => {
  // `filterIndividuals("ind05")` instead of `["ind05"]`: a string spread
  // into an array is its letters, so the call would ask for the individuals
  // `i`, `n`, `d`, `0` and `5`, and the user would read that `i` is not an
  // individual of the variants. The message says what to write instead, as
  // the one of `fields` does for one field.
  const variants = many();

  assert.throws(
    () => variants.filterIndividuals("ind05" as unknown as string[]),
    (error: unknown) =>
      error instanceof Error &&
      error.message.includes("`individuals`") &&
      error.message.includes('["ind05"]'),
  );

  assert.deepEqual(variants.steps, []);
  variants.free();
});
