/**
 * That a range of a file that comes back shorter than the one popnei asked
 * for ends the pass with an error, and does not pass for the end of the
 * file.
 *
 * A browser gives a short range when the file changed on disk after the
 * page got its handle. A reader that took those bytes for the end of the
 * file would give the variants it had read so far and say nothing, which is
 * an analysis of part of a dataset reported as the whole of it. So the case
 * makes the browser do it: `Blob.prototype.slice` is replaced, for the
 * length of one pass, by one that hands back one byte less than the range
 * it was asked for.
 *
 * The file is `tests/reference/vcf/many.vcf`, 117346 bytes, which is
 * smaller than the 4 MiB of a range, so the pass asks for the whole of it
 * in one call and gets 117345 bytes. Those three numbers are what the
 * message has to name: the range popnei asked for, where it asked for it
 * and how much came back.
 *
 * The source is opened before the patch is put on, so the reading that
 * fails is a pass over the variants and not the reading of the header at
 * `openVcf`, which is the case a user would meet: they picked the file,
 * popnei read its header, and the file changed while an analysis ran.
 */

import { openVcf } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { pickedFile } from "../picked_file.ts";

/**
 * The bytes of `many.vcf`, which is the range the pass asks for and where,
 * and the bytes the patched `slice` gives back.
 */
const NUM_BYTES = 117346;
const AT = 0;
const NUM_BYTES_GIVEN = NUM_BYTES - 1;

/**
 * What the message of the error has to say, each piece written as popnei
 * writes it: the range that was asked for and what came back.
 */
const THE_RANGE = `${NUM_BYTES} bytes from ${AT}`;
const WHAT_CAME_BACK = `gave ${NUM_BYTES_GIVEN}`;

/**
 * Runs `body` with every `Blob` of this worker giving one byte less than
 * the range it is asked for, and gives back what `body` threw, or
 * `undefined` when it threw nothing.
 *
 * The patch is taken off again whatever happens, so a case that fails here
 * leaves the next case of this worker a browser that reads whole ranges.
 */
function whileEveryRangeComesBackShort(body: () => void): unknown {
  const wholeRange = Blob.prototype.slice;
  Blob.prototype.slice = function oneByteLess(
    this: Blob,
    start?: number,
    end?: number,
    contentType?: string,
  ): Blob {
    const from = start ?? 0;
    const to = end ?? this.size;
    return wholeRange.call(this, from, Math.max(from, to - 1), contentType);
  };
  try {
    body();
    return undefined;
  } catch (error: unknown) {
    return error;
  } finally {
    Blob.prototype.slice = wholeRange;
  }
}

/**
 * Opens a `File` of `many.vcf`, makes a pass over it with every range one
 * byte short, and asserts the error and its message.
 */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/vcf/many.vcf");
  const file = pickedFile("many.vcf", bytes);
  const variants = openVcf(file, { onlyPassed: false });
  let numVars = 0;
  try {
    const thrown = whileEveryRangeComesBackShort(() => {
      for (const block of variants.iterBlocks()) {
        numVars += block.numVars;
      }
    });
    if (thrown === undefined) {
      throw new Error(
        `a pass over a file whose ranges come back one byte short read ` +
          `${numVars} variants and threw nothing`,
      );
    }
    if (!(thrown instanceof Error)) {
      throw new Error(
        `the pass threw ${JSON.stringify(thrown)}, which is no \`Error\``,
      );
    }
    for (const said of [THE_RANGE, WHAT_CAME_BACK]) {
      if (!thrown.message.includes(said)) {
        throw new Error(
          `the message of the short range does not say \`${said}\`: ` +
            thrown.message,
        );
      }
    }
  } finally {
    variants.free();
  }
}
