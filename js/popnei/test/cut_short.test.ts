/**
 * The two files popnei refuses and pyNei reads, through `openVcf`: a
 * quality that is not finite, and a VCF that bgzip wrote and that does not
 * end with the mark of its end.
 *
 * `docs/specs/io_vcf.md` has both under "The cases a reader of the rules
 * would not guess". A quality of `nan` is an error because NaN is what a
 * block holds for a variant with no quality; it is one only when the
 * quality is asked for, since a column that is not parsed is not checked.
 * A file that bgzip wrote ends with an empty gzip block of 28 bytes, and
 * without it the file was cut short: the variants that are there are given
 * first and the error comes where the iteration would have ended.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVcf } from "popnei";

import { referenceVcf, vcfOf } from "./reference.ts";

await init();

/** The bytes of `many.vcf.gz`, which bgzip wrote. */
const MANY_GZ = await referenceVcf("many.vcf.gz");

/** The empty block of 28 bytes that bgzip writes at the end of a file. */
const MARK_OF_THE_END = 28;

test("a bgzipped VCF without the mark of its end is refused after its variants", () => {
  const cut = MANY_GZ.slice(0, MANY_GZ.length - MARK_OF_THE_END);
  const variants = openVcf(cut, { onlyPassed: false });
  try {
    let read = 0;
    assert.throws(() => {
      for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
        read += block.numVars;
      }
    }, /cut short/);
    // Every variant of the file is before the 28 bytes that are missing.
    assert.equal(read, 500);
  } finally {
    variants.free();
  }
});

test("a bgzipped VCF cut inside a member is refused after its variants", () => {
  // The members of `many.vcf.gz` end at the bytes 310, 12336, 21876 and
  // 21904, so this cut is inside the third one: the decoder runs out of
  // bytes there, and what a user is told is that the file is cut short.
  const cut = MANY_GZ.slice(0, 21000);
  const variants = openVcf(cut, { onlyPassed: false });
  try {
    let read = 0;
    assert.throws(() => {
      for (const block of variants.iterBlocks({ numVarsPerBlock: 100 })) {
        read += block.numVars;
      }
    }, /cut short/);
    // 480 data lines arrived whole, of which the four blocks of 100 that
    // the pass filled are given: the 80 that were left over are lost with
    // the error, as `docs/specs/block.md` says of `reblock`.
    assert.equal(read, 400);
  } finally {
    variants.free();
  }
});

test("a quality that is not finite is refused, and only when it is asked for", () => {
  const bytes = vcfOf(["chr1\t10\t.\tA\tT\tnan\tPASS\t.\tGT\t0/0\t0/1\t1/1"]);
  const variants = openVcf(bytes);
  try {
    assert.throws(() => {
      for (const block of variants.iterBlocks({ fields: ["qual"] })) {
        block.qual;
      }
    }, /QUAL/);

    const numVars = [];
    for (const block of variants.iterBlocks()) {
      numVars.push(block.numVars);
    }
    assert.deepEqual(numVars, [1]);
  } finally {
    variants.free();
  }
});
