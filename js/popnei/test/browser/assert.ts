/**
 * What a case of `cases/` asserts with. node's `assert` is of node, and a
 * case runs in a browser.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`, and Playwright the ones whose name ends in `.browser.ts`.
 */

/**
 * That `found` is `expected`, both written out when they are not.
 *
 * The two are compared as JSON, so numbers, strings and arrays of them are
 * what it takes. A typed array is made an array of numbers by its caller,
 * since JSON writes it as an object of its indices.
 *
 * @throws {Error} When the two differ, with `what` in front of both, which
 * is the message the page and then the browser test read.
 */
export function assertEqual(
  what: string,
  found: unknown,
  expected: unknown,
): void {
  const foundAsJson = JSON.stringify(found);
  const expectedAsJson = JSON.stringify(expected);
  if (foundAsJson !== expectedAsJson) {
    throw new Error(`${what}: ${foundAsJson} and not ${expectedAsJson}`);
  }
}

/**
 * The bytes of the file `path` of the repository, from the server the
 * browser tests read it over.
 *
 * @throws {Error} When the server does not answer with the file, which is
 * what a path that is not one of the repository gives.
 */
export async function bytesOf(path: string): Promise<Uint8Array> {
  const answer = await fetch(path);
  if (!answer.ok) {
    throw new Error(`${path}: the server answered ${answer.status}`);
  }
  return new Uint8Array(await answer.arrayBuffer());
}
