/**
 * That the entry point everything that is not node gets, `dist/web.js`,
 * loads the WebAssembly and answers.
 *
 * node resolves the name of the package to `dist/node.js`, so this file
 * imports the other entry point by its path. That one leaves the loading
 * to the JavaScript wasm-bindgen generated, which fetches the wasm file
 * from the address it sits at; the `fetch` of node does not open a `file:`
 * address, so the test puts in its place one that reads the file and
 * answers as a server would. What is tested is the entry point, not the
 * fetch.
 */

import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

import { versionOfTheCoreCrate } from "./manifest.ts";

import { init, version } from "../dist/web.js";

globalThis.fetch = (async (address: URL | string) => {
  const wasm = await readFile(new URL(address));
  return new Response(wasm, {
    headers: { "content-type": "application/wasm" },
  });
}) as typeof fetch;

test("the entry point of a page fetches the WebAssembly and answers", async () => {
  await init();
  assert.equal(version(), await versionOfTheCoreCrate());
});
