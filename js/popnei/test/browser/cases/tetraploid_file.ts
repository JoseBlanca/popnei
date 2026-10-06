/**
 * That a VCF the user picked in the page, opened with no ploidy, is read
 * with the ploidy of its first genotype with alleles.
 *
 * With no ploidy, `openVcf` of a `File` reads the file twice before it
 * returns: once for the ploidy and once for the individuals, each an
 * opening pass of its own over the ranges of the file. This is the one test
 * of the second of them, which no test under node can make, since a `File`
 * is read only inside a web worker.
 *
 * `tests/reference/dists/tetraploid.vcf.gz` is the tetraploid dataset of
 * `docs/specs/dists.md`, 200 variants of 12 individuals, and a ploidy of 4
 * is one that no default gives. The numbers are those the tests under node
 * assert of the same file.
 */

import { openVcf } from "../../../dist/web.js";
import { assertEqual, bytesOf } from "../assert.ts";
import { pickedFile } from "../picked_file.ts";

/**
 * Opens a `File` of `tetraploid.vcf.gz` with no ploidy and asserts its
 * ploidy and its variants.
 */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/dists/tetraploid.vcf.gz");
  const file = pickedFile("tetraploid.vcf.gz", bytes);
  const variants = openVcf(file);
  try {
    assertEqual("the ploidy of tetraploid.vcf.gz", variants.ploidy, 4);
    assertEqual(
      "the individuals of tetraploid.vcf.gz",
      variants.numIndividuals,
      12,
    );
    let numVars = 0;
    let numAlleles = 0;
    for (const block of variants.iterBlocks()) {
      numVars += block.numVars;
      numAlleles += block.gts.length;
    }
    assertEqual("the variants of tetraploid.vcf.gz", numVars, 200);
    // Four alleles for each of the 12 individuals of each of the 200
    // variants.
    assertEqual("the alleles of tetraploid.vcf.gz", numAlleles, 9600);
  } finally {
    variants.free();
  }
}
