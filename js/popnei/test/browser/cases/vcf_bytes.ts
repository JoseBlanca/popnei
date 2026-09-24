/**
 * That popnei reads a VCF held as an array of bytes inside a web worker
 * and gives there the block the tests under node give.
 *
 * The file is `tests/reference/vcf/many.vcf`, 500 variants of 50
 * individuals, and the numbers are the literals of the tests under node
 * over that same file: the block of 100 variants of `test/vcf.test.ts`,
 * and the individuals, the first five positions and the genotypes of
 * `ind00` of `test/filter_individuals.test.ts`, which reads them from the
 * table of `docs/specs/io_vcf.md`. `test/vcf.test.ts` asserts no column of
 * that block, so the chromosomes, the positions and the genotypes here are
 * the ones of the second file.
 *
 * What it says about the browser is that the WebAssembly loads in a module
 * worker and answers with the numbers the same code answers with under
 * node. The reading of a file by ranges is another case, of work package 3
 * of `docs/plans/js-sources.md`.
 */

import { openVcf } from "../../../dist/web.js";
import { assertEqual, bytesOf } from "../assert.ts";

/** The 50 individuals of `many.vcf`, and its first block of 100 variants. */
const NUM_INDIVIDUALS = 50;
const NUM_VARS_OF_THE_FIRST_BLOCK = 100;

/** The chromosome and the position of its first five variants. */
const FIRST_CHROMS = ["chr1", "chr1", "chr1", "chr1", "chr1"];
const FIRST_POSITIONS = [1000, 1037, 1074, 1111, 1148];

/**
 * The genotypes of `ind00` at the first of those variants and at the third,
 * which the VCF writes as `1/1` and `2|1`.
 */
const GTS_OF_IND00_AT_1000 = [1, 1];
const GTS_OF_IND00_AT_1074 = [2, 1];

/** Reads `many.vcf` from its bytes and asserts its first block. */
export async function run(): Promise<void> {
  const variants = openVcf(await bytesOf("/tests/reference/vcf/many.vcf"));
  try {
    assertEqual("the individuals", variants.numIndividuals, NUM_INDIVIDUALS);
    for (const block of variants.iterBlocks({
      numVarsPerBlock: NUM_VARS_OF_THE_FIRST_BLOCK,
    })) {
      assertEqual(
        "the variants of the first block",
        block.numVars,
        NUM_VARS_OF_THE_FIRST_BLOCK,
      );
      assertEqual(
        "the chromosomes of its first five variants",
        block.chrom?.slice(0, FIRST_CHROMS.length),
        FIRST_CHROMS,
      );
      assertEqual(
        "the positions of its first five variants",
        [...(block.pos ?? []).slice(0, FIRST_POSITIONS.length)],
        FIRST_POSITIONS,
      );
      // The genotypes of a block are the variants one after another, and
      // each variant its individuals one after another, `ind00` first.
      const ofAVariant = block.numIndividuals * block.ploidy;
      assertEqual(
        "the genotypes of ind00 at chr1 1000",
        [...block.gts.subarray(0, block.ploidy)],
        GTS_OF_IND00_AT_1000,
      );
      assertEqual(
        "the genotypes of ind00 at chr1 1074",
        [...block.gts.subarray(2 * ofAVariant, 2 * ofAVariant + block.ploidy)],
        GTS_OF_IND00_AT_1074,
      );
      // The first block is the one this case is about, and leaving the
      // iteration here gives back the memory of wasm of its pass.
      break;
    }
  } finally {
    variants.free();
  }
}
