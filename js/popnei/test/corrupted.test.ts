/**
 * A bgzipped VCF whose bytes were damaged after bgzip wrote it.
 *
 * A review of the VCF reader found one change of two bytes of
 * `tests/reference/vcf/many.vcf.gz` that made the whole file read as no
 * variant and no error: the length of the extra field of the second member,
 * `06 00`, read as `44 54`, which says that the member is longer than it is.
 * The owner decided on 21 September 2026 that such a file is an error
 * however improbable the damage, so the reader reads a bgzipped file by the
 * size that each member states and checks what comes out of it.
 *
 * JavaScript has one exception for everything a library refuses, so this is
 * an `Error` with the message the core gives it, where Python has the
 * `OSError` of a file that cannot be read. `tests/test_corrupted_bgzip.py`
 * reads the same file.
 */

import assert from "node:assert/strict";
import { test } from "node:test";

import { init, openVcf } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/**
 * Where the length of the extra field of the second member is, and what
 * bgzip wrote there: 6 bytes, the one field of a bgzip member. The member
 * starts at the byte 310 of the file.
 */
const THE_FIRST_BYTE = 320;
const AS_BGZIP_WROTE_THEM = [0x06, 0x00];
const AS_THE_REVIEW_CHANGED_THEM = [0x44, 0x54];

/**
 * How many variants are in the members before the damaged one: the first
 * member of `many.vcf.gz` holds the header and no whole data line.
 */
const VARIANTS_BEFORE_THE_DAMAGED_MEMBER = 0;

/** `many.vcf.gz` with those two bytes changed. */
async function withTheDamagedMember(): Promise<Uint8Array> {
  const bytes = await referenceVcf("many.vcf.gz");
  assert.deepEqual(
    Array.from(bytes.slice(THE_FIRST_BYTE, THE_FIRST_BYTE + 2)),
    AS_BGZIP_WROTE_THEM,
  );
  bytes.set(AS_THE_REVIEW_CHANGED_THEM, THE_FIRST_BYTE);
  return bytes;
}

test("a bgzipped VCF whose member is damaged throws and gives no empty file", async () => {
  const variants = openVcf(await withTheDamagedMember(), {
    onlyPassed: false,
  });
  try {
    let read = 0;
    assert.throws(
      () => {
        for (const block of variants.iterBlocks()) {
          read += block.numVars;
        }
      },
      // A user who gets this looks at the member with `xxd -s 310`, and
      // then fetches the file again.
      /member 2, which starts at the byte 310 .* is corrupted/s,
    );
    assert.equal(read, VARIANTS_BEFORE_THE_DAMAGED_MEMBER);
  } finally {
    variants.free();
  }
});
