/**
 * `writeVcf`: the variants of a pass written as the bytes of a VCF.
 *
 * The test is the one `docs/specs/io_vcf.md` gives to TypeScript under
 * "How it is verified" of the writer: `many.vcf` read from a `Uint8Array`
 * with `onlyPassed` false and written with the default, bgzipped, whose
 * bytes decompress to those of the file, and written with `{bgzip: false}`,
 * which gives them as they are. The members are decompressed by node's
 * zlib, which knows nothing of popnei.
 */

import assert from "node:assert/strict";
import { test } from "node:test";
import { gunzipSync } from "node:zlib";

import { init, openVcf, writeVcf } from "popnei";

import { referenceVcf } from "./reference.ts";

await init();

/** The empty member of 28 bytes that ends a file bgzip wrote. */
const THE_EMPTY_MEMBER = Uint8Array.from(
  Buffer.from("1f8b08040000000000ff0600424302001b0003000000000000000000", "hex"),
);

test("writeVcf of many.vcf read with every line gives its bytes, bgzipped by default and plain", async () => {
  const many = await referenceVcf("many.vcf");
  const variants = openVcf(many, { onlyPassed: false });
  try {
    const bgzipped = writeVcf(variants);
    assert.equal(bgzipped.passStats.numVars, 500);
    assert.deepEqual(new Uint8Array(gunzipSync(bgzipped.bytes)), many);
    assert.deepEqual(bgzipped.bytes.subarray(-28), THE_EMPTY_MEMBER);

    const plain = writeVcf(variants, { bgzip: false });
    assert.equal(plain.passStats.numVars, 500);
    assert.deepEqual(plain.bytes, many);
  } finally {
    variants.free();
  }
});

test("writeVcf refuses a bgzip that is not a boolean", async () => {
  const variants = openVcf(await referenceVcf("many.vcf"));
  try {
    assert.throws(
      // @ts-expect-error: the option is a boolean, and a page may give it
      // what it read from a form.
      () => writeVcf(variants, { bgzip: "yes" }),
      /`bgzip` is true or false/,
    );
  } finally {
    variants.free();
  }
});
