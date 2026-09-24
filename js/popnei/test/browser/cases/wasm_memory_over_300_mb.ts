/**
 * That a file of a few hundred MB is read without ever being in the memory
 * of wasm whole: one pass over a `File` of 299994147 bytes gives the
 * 1285000 variants of it and leaves the memory of the module under
 * 16777216 bytes, 16 MiB.
 *
 * The file is 17.9 times that bound, so a build that took it into the
 * memory of wasm whole, which is what an application does today when it
 * hands `openVcf` an array of bytes, cannot pass this. That was run and not
 * only reasoned: on 24 September 2026 `NUM_BYTES_PER_RANGE` of
 * `crates/popnei-js/src/source.rs` was set to 536870912 bytes, more than
 * the file, so the pass read it as one range, and this case failed with the
 * memory of wasm at 601423872 bytes, 35.8 times the bound. The constant was
 * put back at 4194304 and the wasm built again.
 *
 * What a pass leaves in that memory is 14155776 bytes at a range of
 * 4194304, the same byte count in each of the four rounds of the
 * measurement of `docs/reports/js-sources-measurement.md` of 24 September
 * 2026. It is three range-sized blocks and 1572864 bytes that do not grow
 * with the range: the pass holds the range it is reading and the one before
 * it, which is freed once the new one is in its place, and the third block
 * is the range that the reading of the header allocated at `openVcf` and
 * freed there.
 *
 * So the bound is 1.19 times what was measured, which is close enough to
 * catch a change that keeps one range more alive, 18350080 bytes, and a
 * range of 6291456 bytes, 20447232. The second of those two was run on
 * 24 September 2026: with `NUM_BYTES_PER_RANGE` at 6291456 this case failed
 * with the memory of wasm at 20447232 bytes. The first is arithmetic over
 * the three blocks above and was not run.
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
import { assertEqual } from "../assert.ts";
import { fileOfManyVcfRepeated } from "../picked_file.ts";

/** How many times the body of `many.vcf` goes after its header. */
const NUM_COPIES = 2570;

/**
 * How many bytes the file holds and how many variants it gives: 617 plus
 * 2570 times 116729, and 2570 times the 500 variants of `many.vcf`, the 25
 * that failed their FILTER among them, which `onlyPassed` false asks for.
 */
const NUM_BYTES = 299994147;
const NUM_VARS = 1285000;

/**
 * The most the memory of wasm may hold when the pass has ended, 16 MiB,
 * which is 1.19 times the 14155776 bytes the report named above measured.
 */
const MOST_BYTES_OF_WASM = 16777216;

/**
 * Passes once over a `File` of the body of `many.vcf` repeated and asserts
 * its variants and what the memory of wasm came to.
 *
 * @throws {Error} When the file is not of the size these numbers are of,
 * when the pass gives another number of variants, and when the memory of
 * wasm passed the bound.
 */
export async function run(): Promise<void> {
  const what = "the body of many.vcf 2570 times as a File";
  const file = await fileOfManyVcfRepeated(NUM_COPIES);
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
