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
