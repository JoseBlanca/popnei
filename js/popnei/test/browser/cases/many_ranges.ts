/**
 * That a file of more than one range gives the variants of every range of
 * it, and not the variants of its first range.
 *
 * Every other file of these cases is smaller than the 4 MiB of a range, so
 * a pass over it asks the browser for one range and stops: nothing in them
 * reads a second range, and nothing in them crosses from one range into the
 * next. A `Blob.prototype.slice` that gave the end of the first range as
 * the end of the file passed all of them, which is the silent wrong result
 * that a short range is refused for.
 *
 * The file is the header of `tests/reference/vcf/many.vcf`, 617 bytes, with
 * the 116729 bytes of its body after it 120 times: 14008097 bytes and 60000
 * variants of 50 individuals, read as four ranges of 4194304 bytes and one
 * of 1425185. The pass is told how far it has got at its first read, at
 * each of the three reads that bring the bytes it read to a range, and once
 * more at the end of the run, so its five calls are what says that the five
 * ranges were read.
 */

import { openVcf } from "../../../dist/web.js";
import { assertEqual } from "../assert.ts";
import { fileOfManyVcfRepeated } from "../picked_file.ts";

/** How many times the body of `many.vcf` goes after its header. */
const NUM_COPIES = 120;

/**
 * How many bytes the file holds and how many variants it gives: 617 plus
 * 120 times 116729, and 120 times the 500 variants of `many.vcf`, the 25
 * that failed their FILTER among them, which `onlyPassed` false asks for.
 */
const NUM_BYTES = 14008097;
const NUM_VARS = 60000;

/**
 * What the pass has read at each of the calls it makes: nothing at its
 * first read, a range more at each of the next three, and the whole file at
 * the call that ends the run.
 */
const BYTES_READ = [0, 4194304, 8388608, 12582912, 14008097];

/** How many variants a block of the pass holds. */
const NUM_VARS_PER_BLOCK = 1000;

/**
 * Opens a `File` of the body of `many.vcf` repeated and asserts the
 * variants of every range of it and the bytes each call was told.
 */
export async function run(): Promise<void> {
  const what = "the body of many.vcf 120 times as a File";
  const file = await fileOfManyVcfRepeated(NUM_COPIES);
  assertEqual(`${what}: the bytes of the file`, file.size, NUM_BYTES);

  const variants = openVcf(file, { onlyPassed: false });
  const bytesRead: number[] = [];
  variants.onProgress((progress) => {
    bytesRead.push(progress.bytesRead);
  });
  let numVars = 0;
  try {
    for (const block of variants.iterBlocks({
      numVarsPerBlock: NUM_VARS_PER_BLOCK,
    })) {
      numVars += block.numVars;
    }
  } finally {
    variants.free();
  }
  assertEqual(`${what}: the variants of the pass`, numVars, NUM_VARS);
  assertEqual(`${what}: the bytes read at each call`, bytesRead, BYTES_READ);
}
