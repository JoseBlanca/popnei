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
 * The file is `tests/reference/vcf/many.vcf`, 117346 bytes. What the case
 * asserts is not a size: the patched `slice` writes down the range it was
 * asked for, and the message has to name that range, where it was asked
 * for and one byte less than it as what came back. So the case holds
 * whatever size of range popnei reads by, which the measurement of "Speed"
 * of `docs/specs/js_sources.md` sets.
 *
 * The source is opened before the patch is put on, so the reading that
 * fails is a pass over the variants and not the reading of the header at
 * `openVcf`, which is the case a user would meet: they picked the file,
 * popnei read its header, and the file changed while an analysis ran.
 */

import { openVcf } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { pickedFile } from "../picked_file.ts";

/** The range the patched `slice` was asked for, in the order of the asks. */
interface Asked {
  at: number;
  numBytes: number;
}

/**
 * Runs `body` with every `Blob` of this worker giving one byte less than
 * the range it is asked for, and gives back the ranges it was asked for and
 * what `body` threw, or `undefined` when it threw nothing.
 *
 * The patch is taken off again whatever happens, so a case that fails here
 * leaves the next case of this worker a browser that reads whole ranges.
 */
function whileEveryRangeComesBackShort(body: () => void): {
  asked: Asked[];
  thrown: unknown;
} {
  const asked: Asked[] = [];
  const wholeRange = Blob.prototype.slice;
  Blob.prototype.slice = function oneByteLess(
    this: Blob,
    start?: number,
    end?: number,
    contentType?: string,
  ): Blob {
    const from = start ?? 0;
    const to = end ?? this.size;
    asked.push({ at: from, numBytes: to - from });
    return wholeRange.call(this, from, Math.max(from, to - 1), contentType);
  };
  try {
    body();
    return { asked, thrown: undefined };
  } catch (error: unknown) {
    return { asked, thrown: error };
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
    const { asked, thrown } = whileEveryRangeComesBackShort(() => {
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
    // The first range of the pass is the one that came back short, so it is
    // the one the message names. The pass reads no other.
    const first = asked[0];
    if (first === undefined) {
      throw new Error(
        "the pass asked the file for no range, so nothing came back short",
      );
    }
    const theRange = `${first.numBytes} bytes from ${first.at}`;
    const whatCameBack = `gave ${first.numBytes - 1}`;
    for (const said of [theRange, whatCameBack]) {
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
