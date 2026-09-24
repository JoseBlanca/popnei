/**
 * That a file of a few hundred MB is read without ever being in the memory
 * of wasm whole: one pass over a `File` of 299994147 bytes gives the
 * 1285000 variants of it and leaves the memory of the module under
 * 25165824 bytes, 24 MiB.
 *
 * The file is 11.9 times that bound, so a build that took it into the
 * memory of wasm whole, which is what an application does today when it
 * hands `openVcf` an array of bytes, cannot pass this. That was run and not
 * only reasoned: on 24 September 2026 `NUM_BYTES_PER_RANGE` of
 * `crates/popnei-js/src/source.rs` was set to 536870912 bytes, more than
 * the file, so the pass read it as one range, and this case failed with the
 * memory of wasm at 601423872 bytes, 23.9 times the bound. The constant was
 * put back at 4194304 and the wasm built again.
 *
 * Where the two numbers come from: `docs/reports/js-sources-measurement.md`
 * of 24 September 2026 timed one pass over this same file in Chromium and
 * read 14155776 bytes of the memory of wasm at the end of it, the same byte
 * count in each of its four rounds. The bound is 1.78 times that, which
 * leaves room for a consumer that holds more than `calcPerIndividualStats`
 * does. A pass holds about 3.4 times its range of 4194304 bytes, and not
 * one range, because `RangesOfAFile::reads_the_range` builds the next range
 * while the pass still holds the one before it.
 *
 * The memory of a WebAssembly module never shrinks, so what
 * `memory.buffer.byteLength` holds after the pass is the most it held while
 * the pass ran, and one reading at the end is what says it never passed the
 * bound. It is read as `js/popnei/test/vars_memory.test.ts` reads it, from
 * the loader of `wasm/`, which the package of `dist/` loaded before this
 * case was imported. Each case of this directory gets a page and a worker
 * of its own, so nothing else has grown that memory first.
 *
 * The consumer is `calcPerIndividualStats`: it makes one pass, reads the
 * genotypes and adds two counts per individual, so what the memory holds is
 * the reading of the file and not a calculation over it.
 */

import { calcPerIndividualStats, openVcf } from "../../../dist/web.js";
import loadTheWasm from "../../../wasm/popnei.js";
import { assertEqual, bytesOf } from "../assert.ts";
import { aBlockOfItsOwn, pickedFileOfPieces } from "../picked_file.ts";

/** Where the header of `many.vcf` ends, which is the line of its columns. */
const NUM_BYTES_OF_THE_HEADER = 617;

/** How many times the body of the file goes after that header. */
const NUM_COPIES = 2570;

/**
 * How many bytes the file holds and how many variants it gives: 617 plus
 * 2570 times 116729, and 2570 times the 500 variants of `many.vcf`, the 25
 * that failed their FILTER among them, which `onlyPassed` false asks for.
 */
const NUM_BYTES = 299994147;
const NUM_VARS = 1285000;

/**
 * The most the memory of wasm may hold when the pass has ended, 24 MiB, of
 * the report named above.
 */
const MOST_BYTES_OF_WASM = 25165824;

/**
 * Passes once over a `File` of the body of `many.vcf` repeated and asserts
 * its variants and what the memory of wasm came to.
 *
 * @throws {Error} When the file is not of the size these numbers are of,
 * when the pass gives another number of variants, and when the memory of
 * wasm passed the bound.
 */
export async function run(): Promise<void> {
  const many = await bytesOf("/tests/reference/vcf/many.vcf");
  const header = aBlockOfItsOwn(many.subarray(0, NUM_BYTES_OF_THE_HEADER));
  const body = aBlockOfItsOwn(many.subarray(NUM_BYTES_OF_THE_HEADER));
  const pieces = [header];
  for (let copy = 0; copy < NUM_COPIES; copy += 1) {
    // The same block of memory each time: the file reads it 2570 times and
    // the worker holds it once, so the 300 MB are the browser's and never
    // the worker's.
    pieces.push(body);
  }
  const what = "the body of many.vcf 2570 times as a File";
  const file = pickedFileOfPieces("many_repeated.vcf", pieces);
  assertEqual(`${what}: the bytes of the file`, file.size, NUM_BYTES);

  const wasm = await loadTheWasm();
  const variants = openVcf(file, { onlyPassed: false });
  let numVars: number;
  try {
    numVars = calcPerIndividualStats(variants).passStats.numVars;
  } finally {
    variants.free();
  }
  const memoryOfWasm = wasm.memory.buffer.byteLength;

  assertEqual(`${what}: the variants of the pass`, numVars, NUM_VARS);
  if (memoryOfWasm > MOST_BYTES_OF_WASM) {
    throw new Error(
      `${what}: the pass left the memory of wasm at ${memoryOfWasm} bytes, ` +
        `more than the ${MOST_BYTES_OF_WASM} it is held to, over a file of ` +
        `${NUM_BYTES} bytes`,
    );
  }
}
