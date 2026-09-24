/**
 * The `File` a case of `cases/` hands popnei, built from bytes it read over
 * http or wrote itself.
 *
 * A `File` is the handle a page gets when the user picks a file in a form
 * or drops one on it: it carries the name and the size of the file and
 * gives any range of its bytes. A browser test cannot make the user pick
 * anything, so it builds the handle from bytes, and what popnei then does
 * with it, asking for one range at a time, is the same.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`, and Playwright the ones whose name ends in `.browser.ts`.
 */

import { bytesOf } from "./assert.ts";

/** Where the header of `tests/reference/vcf/many.vcf` ends, which is the
 * line of its columns. */
const NUM_BYTES_OF_THE_HEADER = 617;

/**
 * The `File` called `name` holding `bytes`.
 *
 * The bytes are copied into a block of memory of their own, because an
 * array of a worker that shares its memory with the page sits on a
 * `SharedArrayBuffer`, which a `File` does not take.
 */
export function pickedFile(name: string, bytes: Uint8Array): File {
  return pickedFileOfPieces(name, [aBlockOfItsOwn(bytes)]);
}

/**
 * The `File` called `name` holding the pieces of `pieces`, one after
 * another.
 *
 * It is how a case builds a file of several ranges without holding its
 * bytes as many times as it repeats them: one piece given twice is one
 * block of memory that the `File` reads twice.
 */
export function pickedFileOfPieces(name: string, pieces: ArrayBuffer[]): File {
  return new File(pieces, name);
}

/**
 * `bytes` copied into a block of memory of their own, which is what a
 * `File` takes and what the bytes of a `Uint8Array` of TypeScript are not
 * known to be.
 */
export function aBlockOfItsOwn(bytes: Uint8Array): ArrayBuffer {
  const buffer = new ArrayBuffer(bytes.byteLength);
  new Uint8Array(buffer).set(bytes);
  return buffer;
}

/**
 * The `File` of the header of `tests/reference/vcf/many.vcf` with the body
 * of that file after it `numCopies` times: 617 + 116729 * `numCopies`
 * bytes, and 500 * `numCopies` variants of 50 diploid individuals, the 25
 * of each copy that failed their FILTER among them.
 *
 * The body goes in as one block of memory that the `File` reads again for
 * each copy, so a file of 300 MB costs the worker the 117346 bytes of
 * `many.vcf` and the browser the rest. A case that opens one asserts the
 * size of the file and the variants of the pass against its own literals,
 * which is what says the file is the one it asked for.
 */
export async function fileOfManyVcfRepeated(numCopies: number): Promise<File> {
  const many = await bytesOf("/tests/reference/vcf/many.vcf");
  const header = aBlockOfItsOwn(many.subarray(0, NUM_BYTES_OF_THE_HEADER));
  const body = aBlockOfItsOwn(many.subarray(NUM_BYTES_OF_THE_HEADER));
  const pieces = [header];
  for (let copy = 0; copy < numCopies; copy += 1) {
    pieces.push(body);
  }
  return pickedFileOfPieces("many_repeated.vcf", pieces);
}
