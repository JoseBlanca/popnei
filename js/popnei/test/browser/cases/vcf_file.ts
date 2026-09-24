/**
 * That popnei reads a VCF the user picked in the page, one range of bytes
 * at a time, and gives the variants it gives when the whole file is in the
 * memory of wasm.
 *
 * The `File` is built in the worker from the bytes of
 * `tests/reference/vcf/many.vcf`, which is the handle a page gets from a
 * form or a drop with the bytes of the machine behind it. What popnei then
 * does with it is what no test under node can run: it asks the browser for
 * the ranges it needs through `FileReaderSync`, which exists only inside a
 * web worker.
 *
 * The numbers are those of `../many.ts`, the literals of the tests under
 * node. The file is 117346 bytes and a range is 4 MiB, so this pass reads
 * one range and no more: what it says is that a range arrives in the memory
 * of wasm as the bytes of the file, not that the reader crosses from one
 * range into the next.
 */

import { openVcf } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { assertsTheVariantsOfMany } from "../many.ts";
import { pickedFile } from "../picked_file.ts";

/** Opens a `File` of `many.vcf` and asserts the variants of the file. */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/vcf/many.vcf");
  const file = pickedFile("many.vcf", bytes);
  const variants = openVcf(file, { onlyPassed: false });
  try {
    assertsTheVariantsOfMany("many.vcf as a File", variants);
  } finally {
    variants.free();
  }
}
