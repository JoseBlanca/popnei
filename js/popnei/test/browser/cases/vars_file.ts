/**
 * That a vars file the user picked in the page gives the variants that were
 * written into it, which is the check that popnei can seek inside a file of
 * the page.
 *
 * The case writes the vars file itself: it opens `many.vcf` from its bytes,
 * calls `writeVars`, which gives the whole file as an array of bytes
 * because a tab has no filesystem, and makes a `File` of those bytes, which
 * is what a user who saved that download and picked it again would hand the
 * page.
 *
 * A vars file is one arrow IPC file, and its reader goes to the footer at
 * the end of the file, reads where each batch of variants sits and then
 * goes back to the batches, so reading one means moving where the pass is
 * without reading forward from there. The reader of a VCF never does, so
 * this case is the only one that says the `Seek` over the ranges of a file
 * lands where it was asked to.
 *
 * The numbers are those of `../many.ts`, the literals of the tests under
 * node: `writeVars` writes every variant its source gives, so the file
 * holds the same 500 variants of the same 50 individuals.
 */

import { openVars, openVcf, writeVars } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { assertsTheVariantsOfMany } from "../many.ts";
import { pickedFile } from "../picked_file.ts";

/**
 * The bytes of a vars file of every variant of the VCF in `bytes`, the ones
 * that failed their FILTER among them.
 *
 * @throws {Error} When the bytes are not a VCF popnei can read.
 */
function varsFileOf(bytes: Uint8Array): Uint8Array {
  const ofTheVcf = openVcf(bytes, { onlyPassed: false });
  try {
    return writeVars(ofTheVcf).bytes;
  } finally {
    ofTheVcf.free();
  }
}

/**
 * Writes a vars file of `many.vcf`, opens a `File` of it and asserts the
 * variants of the file.
 */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/vcf/many.vcf");
  const file = pickedFile("many.vars", varsFileOf(bytes));
  const variants = openVars(file);
  try {
    assertsTheVariantsOfMany("a vars file of many.vcf as a File", variants);
  } finally {
    variants.free();
  }
}
