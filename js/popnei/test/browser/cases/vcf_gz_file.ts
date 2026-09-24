/**
 * That a gzipped VCF the user picked in the page gives the variants of the
 * same VCF plain, which is the ranges of a file read through the
 * decompressor.
 *
 * `tests/reference/vcf/many.vcf.gz` is `many.vcf` as bgzip wrote it, 21904
 * bytes of four gzip members against the 117346 of the plain file, and the
 * variants inside are the same 500. The decompressor reads its input
 * through the same `Read` that the plain file goes through, so what this
 * adds to `vcf_file.ts` is that the bytes of a range are handed to it and
 * not to the parser of the lines.
 *
 * The numbers are those of `../many.ts`, the literals of the tests under
 * node.
 */

import { openVcf } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { assertsTheVariantsOfMany } from "../many.ts";
import { pickedFile } from "../picked_file.ts";

/** Opens a `File` of `many.vcf.gz` and asserts the variants of the file. */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/vcf/many.vcf.gz");
  const file = pickedFile("many.vcf.gz", bytes);
  const variants = openVcf(file, { onlyPassed: false });
  try {
    assertsTheVariantsOfMany("many.vcf.gz as a File", variants);
  } finally {
    variants.free();
  }
}
