/**
 * The server the browser tests read the repository from: the built package
 * of `js/popnei/dist/`, the WebAssembly of `js/popnei/wasm/`, the page and
 * the worker of this directory, and the reference files of
 * `tests/reference/`.
 *
 * A module worker and the WebAssembly are loaded with `fetch`, which does
 * not open a `file:` address, so the files a browser test reads are served
 * over http from the root of the repository, at the port the first
 * argument of the command gives. `playwright.config.ts` starts it and
 * stops it.
 *
 * A file whose name ends in `.ts` is served as JavaScript with its types
 * stripped, which is what node does with the tests of `test/`, so the page,
 * the worker and the cases are TypeScript and the browser is given
 * JavaScript. Nothing else of the file is changed: the addresses it imports
 * from are the ones it holds.
 *
 * It is not a test file: node's test runner runs the files whose name ends
 * in `.test.ts`, and Playwright the ones whose name ends in `.browser.ts`.
 */

import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { stripTypeScriptTypes } from "node:module";
import { fileURLToPath } from "node:url";

/** The root of the repository, four directories above this file. */
const REPOSITORY = new URL("../../../../", import.meta.url);

/** What a browser is told each kind of file is. */
const CONTENT_TYPES: Record<string, string> = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".ts": "text/javascript; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".wasm": "application/wasm",
  ".vcf": "text/plain; charset=utf-8",
  ".gz": "application/gzip",
  ".vars": "application/octet-stream",
};

/**
 * The bytes to answer the address `path` with, and what they are.
 *
 * @throws {Error} When the path leaves the repository, which is what a
 * page that asked for a file of the machine outside it would look like.
 */
async function fileOf(path: string): Promise<{ body: Uint8Array; type: string }> {
  const file = new URL(path.replace(/^\/+/, ""), REPOSITORY);
  if (!file.href.startsWith(REPOSITORY.href)) {
    throw new Error(`the address \`${path}\` leaves the repository`);
  }
  const extension = /\.[^./]+$/.exec(file.pathname)?.[0] ?? "";
  const type = CONTENT_TYPES[extension] ?? "application/octet-stream";
  const bytes = await readFile(fileURLToPath(file));
  if (extension !== ".ts") {
    return { body: new Uint8Array(bytes), type };
  }
  const javascript = stripTypeScriptTypes(bytes.toString("utf8"));
  return { body: new TextEncoder().encode(javascript), type };
}

const port = Number(process.argv[2]);
if (!Number.isInteger(port) || port <= 0) {
  throw new Error(`the port of the server is \`${process.argv[2]}\``);
}

createServer((request, response) => {
  const path = new URL(request.url ?? "/", "http://127.0.0.1").pathname;
  fileOf(path).then(
    ({ body, type }) => {
      response.writeHead(200, { "content-type": type });
      response.end(body);
    },
    (error: unknown) => {
      // A browser that asked for a file that is not there gets the message
      // and the test that reads it says which address failed.
      response.writeHead(404, { "content-type": "text/plain; charset=utf-8" });
      response.end(`${path}: ${error instanceof Error ? error.message : String(error)}`);
    },
  );
}).listen(port, "127.0.0.1");
