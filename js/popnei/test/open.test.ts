/**
 * That opening a VCF reads its header and asks for no block.
 *
 * A block of a VCF holds as many variants as the individuals of the file
 * make room for, and the genotypes of one are its variants times its
 * individuals times its ploidy. In wasm a count of things is 32 bits and
 * holds 4295 million, so a file of 170000 individuals read with the ploidy
 * 255 has no block of the size popnei chooses, 100 variants: those are 4335
 * million genotypes. `openVcf` reads the header and nothing else, so such a
 * file is opened, its individuals are read, and the size of its blocks is
 * the user's to choose.
 *
 * It is the case that a review of the binding crates found refused at
 * `openVcf`, and `docs/specs/io_vcf.md` has the rule under "The Rust
 * interface".
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVcf } from "popnei";

await init();

/** How many individuals the header of this test names. */
const MANY_INDIVIDUALS = 170000;

/** The ploidy it is read with, the largest a reader of popnei takes. */
const LARGEST_PLOIDY = 255;

/** A VCF of that many individuals and no variant. */
function headerOfManyIndividuals(numIndividuals: number): Uint8Array {
  const names = Array.from(
    { length: numIndividuals },
    (_unused, individual) => `ind${individual}`,
  );
  const columns = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    ...names,
  ];
  return new TextEncoder().encode(
    ["##fileformat=VCFv4.4", columns.join("\t"), ""].join("\n"),
  );
}

test("a header whose blocks of the size popnei chooses would not fit is opened", () => {
  const bytes = headerOfManyIndividuals(MANY_INDIVIDUALS);
  const variants = openVcf(bytes, { ploidy: LARGEST_PLOIDY });
  try {
    assert.equal(variants.numIndividuals, MANY_INDIVIDUALS);
    assert.equal(variants.ploidy, LARGEST_PLOIDY);
    assert.equal(variants.individuals[0], "ind0");
  } finally {
    variants.free();
  }
});

test("the blocks of such a file are read in a size that fits", async () => {
  const bytes = headerOfManyIndividuals(MANY_INDIVIDUALS);
  const variants = openVcf(bytes, { ploidy: LARGEST_PLOIDY });
  try {
    // The file has a header and no variant, so what this asks of the reader
    // is one block of 10 variants, 433 million genotypes, and the reader
    // gives no block because there is no variant to put in one.
    const blocks: number[] = [];
    for (const block of variants.iterBlocks({ numVarsPerBlock: 10 })) {
      blocks.push(block.numVars);
    }
    assert.deepEqual(blocks, []);

    // A size that does not fit is refused, which is what opening the file
    // does not do: 100 variants of this file are 4335 million genotypes.
    const refused: number[] = [];
    assert.throws(() => {
      for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
        refused.push(block.numVars);
      }
    }, /memory/);
    assert.deepEqual(refused, []);
  } finally {
    variants.free();
  }
});
