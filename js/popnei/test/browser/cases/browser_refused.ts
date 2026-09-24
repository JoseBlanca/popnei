/**
 * That what the browser threw when it refused a range reaches the user as
 * the sentence the browser wrote, and not as the form Rust prints a value
 * of JavaScript it cannot read in.
 *
 * A browser refuses a range when the file was moved or changed on disk
 * after the page got its handle, and what `FileReaderSync` throws then is a
 * `DOMException` called `NotReadableError`, whose `message` is the sentence
 * the user is meant to read. An application that throws its own `Error`
 * from a patched `Blob` is the other value popnei can be given here, and
 * its message is the sentence too.
 *
 * The case makes the browser refuse: `Blob.prototype.slice`, for the length
 * of one pass, throws instead of giving the range. The file is
 * `tests/reference/vcf/many.vcf`, and the source is opened before the patch
 * is put on, so what fails is a pass over the variants and not the reading
 * of the header at `openVcf`.
 */

import { openVcf } from "../../../dist/web.js";
import { bytesOf } from "../assert.ts";
import { pickedFile } from "../picked_file.ts";

/** The sentence of each of the two values the patched `slice` throws. */
const WHAT_THE_BROWSER_SAID = "The requested file could not be read.";
const WHAT_AN_APPLICATION_SAID = "the file is not on this disk any more";

/**
 * How Rust writes a value of JavaScript it has no sentence for, which is
 * what the message must not hold.
 */
const THE_FORM_OF_RUST = "JsValue(";

/**
 * Runs `body` with every `Blob` of this worker throwing `thrown` in place
 * of the range it is asked for, and gives back what `body` threw.
 *
 * The patch is taken off again whatever happens, so a case that fails here
 * leaves the next case of this worker a browser that gives its ranges.
 */
function whileEveryRangeIsRefused(thrown: unknown, body: () => void): unknown {
  const wholeRange = Blob.prototype.slice;
  Blob.prototype.slice = function refuses(): Blob {
    throw thrown;
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
 * That a pass over `many.vcf` whose ranges are refused with `thrown` fails
 * with a message that holds `said` and not the form of Rust.
 *
 * @throws {Error} When the pass gave its variants, when what it threw is no
 * `Error`, and when the message does not say what the browser said.
 */
function theMessageOfARangeRefusedWith(
  what: string,
  thrown: unknown,
  said: string,
  bytes: Uint8Array,
): void {
  const variants = openVcf(pickedFile("many.vcf", bytes), {
    onlyPassed: false,
  });
  try {
    const failure = whileEveryRangeIsRefused(thrown, () => {
      for (const _block of variants.iterBlocks()) {
        // The first block is the one that reads the range.
      }
    });
    if (failure === undefined) {
      throw new Error(`${what}: the pass read the file and threw nothing`);
    }
    if (!(failure instanceof Error)) {
      throw new Error(
        `${what}: the pass threw ${JSON.stringify(failure)}, which is no ` +
          `\`Error\``,
      );
    }
    if (!failure.message.includes(said)) {
      throw new Error(
        `${what}: the message does not say \`${said}\`: ${failure.message}`,
      );
    }
    if (failure.message.includes(THE_FORM_OF_RUST)) {
      throw new Error(
        `${what}: the message gives the value as Rust writes it: ` +
          failure.message,
      );
    }
  } finally {
    variants.free();
  }
}

/** Refuses the ranges of a pass with each of the two values and asserts. */
export async function run(): Promise<void> {
  const bytes = await bytesOf("/tests/reference/vcf/many.vcf");
  theMessageOfARangeRefusedWith(
    "a range the browser refused",
    new DOMException(WHAT_THE_BROWSER_SAID, "NotReadableError"),
    WHAT_THE_BROWSER_SAID,
    bytes,
  );
  theMessageOfARangeRefusedWith(
    "a range an application refused",
    new Error(WHAT_AN_APPLICATION_SAID),
    WHAT_AN_APPLICATION_SAID,
    bytes,
  );
}
